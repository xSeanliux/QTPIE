from Bio import Phylo 
from io import StringIO
from typing import List
from collections import Counter

class QuartetSorter: 
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
        self.assign_internal_labels()
        self.all_leaf_names = map(
            lambda x : x.name, 
            self.tree.find_elements(terminal=True)
        )
        self.polytomies: List[Phylo.BaseTree.Clade] = list(filter(
            lambda x : len(x.clades) > 2,
            self.tree.find_clades(),
        ))
        self.precompute_label_map()
        self.polytomy_quartets = { 
            polytomy.name: Counter() 
            for polytomy in self.polytomies 
        }

    def get_parent(self, node):
        path = [None, self.tree.root] + self.tree.get_path(node)
        if not path:
            raise ValueError(f"Could not find parent of {node=}")
        parent = path[-2]
        return parent

    def precompute_label_map(self):
        poly_label_map = { 
            polytomy.name: {} for polytomy in self.polytomies 
        }
        for poly in self.polytomies:
            parent = self.get_parent(poly)

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
        w: int,
        polytomy_quartets: dict[str, Counter],
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
                polytomy_quartets[pn][relabelled_tuple] += w
                return # a quartet can be in at most one polytomy

    def get_polytomy_quartets(self):
        return self.polytomy_quartets




