import argparse
from Bio import Phylo
import sys
from io import StringIO
import treeswift
from __init__ import QuartetPolytomy

import re

def main():
    parser = argparse.ArgumentParser(description="Process guide tree, quartets, and ASTRAL file.")
    
    parser.add_argument(
        '-g', '--guidetree', 
        type=str, 
        required=True, 
        help="Path to the guide tree file."
    )
    
    parser.add_argument(
        '-f', '--format', 
        type=str, 
        required=False, 
        default='newick',
        help="Guide tree format."
    )
    
    parser.add_argument(
        '-q', '--quartets', 
        type=str, 
        required=True, 
        help="Path to the quartets file."
    )
    
    parser.add_argument(
        '-A', '--astralpath', 
        type=str, 
        default='ASTRAL/astral.5.7.8.jar',
        required=False, 
        help="Path to the ASTRAL executable."
    )

    parser.add_argument(
        '-o', '--output', 
        type=str, 
        required=False, 
        default="./tmp",
        help="Path to output folder."
    )
    
    args = parser.parse_args()
    
    print(f"Guide Tree Path: {args.guidetree} in {args.format} format", file = sys.stderr)
    print(f"Quartets Path: {args.quartets}", file = sys.stderr)
    print(f"ASTRAL File Path: {args.astralpath}", file = sys.stderr)

    qp = QuartetPolytomy(
        tree_file_path = args.guidetree,
        format = args.format
    )
    print("Finished initialising QuartetPolytomy object & precomputation.", file = sys.stderr)
    # update all quartets 
    quartet_re = re.compile(r'\(\((\w+),(\w+)\),\((\w+),(\w+)\)\)\;?')
    with open(args.quartets, "r") as qf:
        quartet_lines = qf.readlines()
        for line in quartet_lines:
            if line == "":
                continue
            m = quartet_re.match(line)
            assert m is not None, f"Quartet parsing error: {line}"
            q = m.groups()
            qp.update_quartet(q)
    print("Finished updating quartets.", file = sys.stderr)
    # run ASTRAL and resolve
    qp.run_resolve(
        output_path = args.output,
        astral_path = args.astralpath
    )

    # clean up extra nodes and remove unifurcations 
    # unfortunately I have to extract the Newick string from Bio.Phylo
    # and use the treeswift suppress_unifurcations function
    for clade in qp.tree.find_clades(terminal=False):
        clade.name = None

    writer = Phylo.NewickIO.Writer([qp.tree])
    res_newick = next(iter((writer.to_strings(plain=True))))
    ts_tree = treeswift.read_tree(res_newick, schema = "newick")
    ts_tree.suppress_unifurcations()
    res_newick = ts_tree.newick()
    print(res_newick)


if __name__ == "__main__":
    main()
