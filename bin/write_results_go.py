#!/usr/bin/env python2
"""
Go-aware adapter for PHYLOWGS write_results.py.

phylowgs-go's trees.zip uses JSON entries (not pickle), so the upstream
util2.TreeReader / pwgsresults.result_generator.ResultGenerator pipeline
can't parse it directly. Each JSON entry, however, already carries exactly
the fields ResultGenerator._summarize_pops() would have produced --
{llh, structure, populations, mut_assignments} -- because phylowgs-go's
summarizePops() (main.go) was written to match that function field-for-field.
So instead of reconstructing a live tree object and re-deriving those fields,
this script reads them straight out of the JSON and hands them to the same
ResultMunger + JsonWriter classes write_results.py uses, giving identical
small-node/superclone/multiprimary filtering and identical output schema
(*.summ.json.gz, *.muts.json.gz, *.mutass.zip).

mutlist (SSM/CNV reference table) comes from phylowgs-go's mutlist.json,
which already matches ResultGenerator._list_mutations()'s output shape.
"""
import argparse
import json
import sys
import zipfile

# pwgsresults isn't pip-installed in the phylowgs container -- it's a plain
# package directory alongside write_results.py at /usr/bin/phylowgs, only
# importable when that dir is on sys.path (which it normally is only for
# scripts invoked from that directory).
sys.path.append('/usr/bin/phylowgs')
from pwgsresults.result_munger import ResultMunger
from pwgsresults.json_writer import JsonWriter


def restricted_float(x):
    x = float(x)
    if not (0 <= x <= 1):
        raise argparse.ArgumentTypeError('%s outside range [0, 1]' % x)
    return x


def _intkeys(d):
    return dict((int(k), v) for k, v in d.items())


def load_go_results(trees_zip_path, mutlist_path):
    with open(mutlist_path) as f:
        mutlist = json.load(f)
    mutlist.pop('dataset_name', None)

    summaries = {}
    mutass = {}
    zf = zipfile.ZipFile(trees_zip_path)
    try:
        for name in zf.namelist():
            if not name.startswith('tree_'):
                continue
            # "tree_<idx>_<llh>"
            idx = int(name.split('_', 2)[1])
            entry = json.loads(zf.read(name))
            summaries[idx] = {
                'llh': entry['llh'],
                'structure': _intkeys(entry['structure']),
                'populations': _intkeys(entry['populations']),
            }
            mutass[idx] = _intkeys(entry.get('mut_assignments', {}))
    finally:
        zf.close()

    return summaries, mutlist, mutass


def main():
    parser = argparse.ArgumentParser(
        description='Write JSON result files from phylowgs-go output, matching write_results.py output schema',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--include-ssm-names', dest='include_ssm_names', action='store_true',
                         help='No-op: retained for CLI parity with write_results.py. phylowgs-go\'s '
                              'mutlist.json already includes SSM names if the input carried them.')
    parser.add_argument('--min-ssms', dest='min_ssms', type=float, default=0.01,
                         help='Minimum number or percent of SSMs to retain a subclone')
    parser.add_argument('--include-multiprimary', dest='include_multiprimary', action='store_true',
                         help='Whether to include multiprimary trees in result')
    parser.add_argument('--keep-superclones', dest='keep_superclones', action='store_true',
                         help='Whether to keep superclones, which are often artifacts of tree sampling')
    parser.add_argument('--max-multiprimary', dest='max_multiprimary', type=restricted_float, default=0.8,
                         help='Maximum proportion of trees that may be multiprimary if '
                              '--include-multiprimary=False. Otherwise an exception is raised if the '
                              'proportion of multiprimary trees exceeds this value.')
    parser.add_argument('dataset_name', help='Name identifying dataset')
    parser.add_argument('tree_file', help='phylowgs-go trees.zip')
    parser.add_argument('mutlist_file', help='phylowgs-go mutlist.json')
    parser.add_argument('tree_summary_output', help='Output file for JSON-formatted tree summaries')
    parser.add_argument('mutlist_output', help='Output file for JSON-formatted list of mutations')
    parser.add_argument('mutass_output', help='Output file for JSON-formatted list of SSMs/CNVs per subclone')
    args = parser.parse_args()

    summaries, mutlist, mutass = load_go_results(args.tree_file, args.mutlist_file)
    params = {}

    munger = ResultMunger(summaries, mutlist, mutass)
    summaries, mutass = munger.remove_small_nodes(args.min_ssms)
    if not args.keep_superclones:
        munger.remove_superclones()
    if not args.include_multiprimary:
        munger.remove_multiprimary_trees(args.max_multiprimary)

    writer = JsonWriter(args.dataset_name)
    writer.write_summaries(summaries, params, args.tree_summary_output)
    writer.write_mutlist(mutlist, args.mutlist_output)
    writer.write_mutass(mutass, args.mutass_output)


if __name__ == '__main__':
    main()
