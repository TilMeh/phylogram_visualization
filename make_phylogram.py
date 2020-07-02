from collections import Counter
import pandas as pd
import sys
import subprocess
import xml.etree.ElementTree as ET
from ete3 import Tree, TreeStyle, TextFace
import annotate_family_tree as aft

species_list_table = aft.get_species_list_from_table("pl_table.csv")
species_list_tree = []
print("Loading Tree...")
t = Tree("largetree_final.tre", format=5)
print("Gathering species(leaf) names...")
for leaf in t.iter_leaves():
    species_list_tree.append(leaf.name.replace("_"," "))
print("Constructing lineage lookup...")
lin_lookup = aft.Lineage_Lookup("fullnamelineage.dmp")
print("Reducing lineage lookup to Spermatophyta...")
lin_lookup.reduce_to_spermatophyta()
print("Constructing dict of species in family...")
species_in_family = aft.get_species_in_family(species_list_tree, lin_lookup)
print("Counting species per family...")
counts_per_family = aft.get_counts_per_family(species_list_table, species_in_family)
print("Pruning Tree to family level...")
aft.prune_to_family(t, species_in_family)
print("Attaching counts to Tree...")
aft.attach_counts_to_tree(t, counts_per_family)
#print("Extending leaf branches...")
#aft.extend_leaf_branches(t)

ts = TreeStyle()
ts.mode = "c"
#ts.optimal_scale_level = "full"
ts.draw_guiding_lines = True
#ts.guiding_lines_type = 0
#ts.guiding_lines_color = "black"
ts.show_leaf_name = False
print("Rendering Tree...")
t.render("family_tree.png", w=10000, h=10000, tree_style=ts)
