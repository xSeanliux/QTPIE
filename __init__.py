from Bio import Phylo 
from io import StringIO
from typing import List
import os
from collections import Counter
from pathlib import Path
from copy import deepcopy


def get_parent(tree, node):
    path = [tree.root] + tree.get_path(node)
    if len(path) <= 1:
        raise ValueError(f"Could not find parent of {node=}, are you sure {node=} is a valid non-root node?")
    parent = path[-2]
    return parent

class QuartetPolytomy: 
    def get_new_name(self):
        self._id += 1
        return f"_{self._id}"

    def reset_name(self):
        self._id = 0

    def assign_internal_labels(self):
        for node in self.tree.find_clades(terminal=False):
            node.name = self.get_new_name()

    def __init__(
        self,
        tree_file_path: str,
        format: str = "newick",
    ):
        self._id = 0
        self.tree = Phylo.read(tree_file_path, format=format)
        self.tree.root_at_midpoint
        self.assign_internal_labels()
        self.all_leaf_names = map(
            lambda x : x.name, 
            self.tree.find_elements(terminal=True)
        )
        self.polytomies = list(filter(
            lambda x : len(x.clades) > 2 
            ,
            self.tree.find_clades(),
        ))
        self.precompute_label_map()
        self.polytomy_quartets = { 
            polytomy.name: Counter() 
            for polytomy in self.polytomies 
        }

    def run_astral(
        self,
        output_path: str,
        astral_path: str,
    ):
        self.write_polytomy_quartets(output_path)
        # run ASTRAL 
        for polytomy in self.polytomies:




    def precompute_label_map(self):
        poly_label_map = { 
            polytomy.name: {} for polytomy in self.polytomies 
        }
        for poly in self.polytomies:
            parent = get_parent(self.tree, poly)

            for child in poly.clades: 
                leaves = list(child.find_elements(terminal=True))
                for leaf in leaves: 
                    poly_label_map[poly.name][leaf.name] = child.name
            if parent is not None:
                leaves = [ # the rest of the leaves
                    l for l in self.tree.find_elements(terminal=True)
                    if l not in poly.find_elements()
                ]
                for leaf in leaves: 
                    poly_label_map[poly.name][leaf.name] = parent.name
            print(len(poly_label_map[poly.name]))
        self.poly_label_map = poly_label_map

    def update_quartet(
        self,
        q: tuple[str, str, str, str],
        w: int = 1,
    ):
        # q has format (a, b, c, d), represents ab|cd
        # please check that a,b,c,d are in the leafset 
        # w is an integer weight
        for poly in self.polytomies: 
            pn = poly.name 
            relabelled_tuple = tuple(map(
                lambda u: self.poly_label_map[pn][u],
                q
            ))
            if len(set(relabelled_tuple)) == 4: # four different children
                self.polytomy_quartets[pn][relabelled_tuple] += w
                return # a quartet can be in at most one polytomy

    def get_polytomy_quartets(self):
        return self.polytomy_quartets

    def write_polytomy_quartets(
        self,
        output_path: str,
    ):
        # For each polytomy with name PN, will write a list of quartets to 
        # output_path/PN/quartets.nwk
        FOLDER = Path(output_path)
        for polytomy in self.polytomies:
            polytomy_folder = FOLDER / polytomy.name
            os.makedirs(polytomy_folder, exist_ok=True)
            # what happens if the leafset is not fully determined? for neighbours a,b,c,d add 
            # ab|cd, ac|bd, and ad|bc so that no bias
            # and make sure each leaf is in at least one quartet.
            with open(polytomy_folder / 'quartets.nwk', 'w') as qf, open(polytomy_folder / 'metadata.csv', 'w') as mdf:
                poly_neighbours = list(set(self.poly_label_map[polytomy.name].values()))

                mdf.write(f"POLYTOMY SIZE {len(poly_neighbours)}\n")

                for i in range(len(poly_neighbours) - 3):
                    a, b, c, d = poly_neighbours[i:i+4]
                    qf.write(f'(({a},{b}),({c},{d}));\n')
                    qf.write(f'(({a},{c}),({b},{d}));\n')
                    qf.write(f'(({a},{d}),({c},{d}));\n')

                for (a, b, c, d), w in self.polytomy_quartets[polytomy.name].items():
                    qf.write(f'(({a},{b}),({c},{d}));\n' * w)
        
    
def resolve_one_polytomy(
    tree_pt, # tree with polytomy
    pt_node_name, # node name of the polytomy to resolve
    resolution_tree, # the resolution tree such that the leafset of this tree is exactly the neighbours of pt_node_name in tree_pt
): 
    """ArithmeticError
    Input:
    tree_pt, tree with polytomy
    pt_node_name, node name of the polytomy to resolve
    resolution_tree, the resolution tree such that the leafset of this tree is exactly the neighbours of pt_node_name in tree_pt
    Output: 
    tree_pt with the resolution tree in place of the original polytomy at pt_node_name. 
    Tree WILL have unifurcations, but they are intended and should be fixed AFTER all polytomies are fixed, as they help to represent nodes that are may no longer present.
    Will not chagne resolution_tree, but CHANGES tree_pt IN PLACE
    """
    # assumptions: tree_pt
    # 
    old_resolution_tree = deepcopy(resolution_tree)
    assert len(list(tree_pt.find_clades(target=pt_node_name))) == 1, f"found multiple nodes with name {pt_node_name}"
    pt_node = next(iter(tree_pt.find_clades(target=pt_node_name)))
    parent = get_parent(tree_pt, pt_node)
    assert parent is not None, "parent of polytomy should not be none."
    assert set([parent] + [ clade.name for clade in pt_node.clades ]) ==\
        set([ leaf.name for leaf in resolution_tree.find_elements(terminal=True) ] ),\
        "Resolution leaf set is not the same as the neighbours of the polytomy"
    # stitch the children the original polytomy on the resolution
    for child in pt_node.clades:
        name = child.name 
        assert len(list(resolution_tree.find_clades(target=name, terminal=True))) == 1, f"Found multiple or no nodes in the resolution with name {name}"
        resolution_child = next(iter(resolution_tree.find_clades(target=name, terminal=True)))
        resolution_child.clades = [child]
        resolution_child.name = pt_node.name
        print(resolution_child, resolution_child.clades)

    # replace the parent -> polytomy node with the tree in tree_pt
    resolution_tree.root_with_outgroup([parent.name])
    resolution_tree.prune(parent.name)
    for i, c in enumerate(parent.clades):
        if c.name == pt_node.name:
            parent.clades[i] = resolution_tree.root
            break

    resolution_tree = old_resolution_tree
    return tree_pt


    
            




