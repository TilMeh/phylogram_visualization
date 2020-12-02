#####################
# IMPORT OPERATIONS #
#####################
from collections import Counter
import pandas as pd
from pathlib import Path
import os.path, time
import subprocess
import xml.etree.ElementTree as ET
import argparse
import coloredlogs, logging
from ete3 import Tree, TreeStyle, TextFace, NCBITaxa

###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>, '\
			 'Tilman Mehl <tilmanmehl@zedat.fu-berlin.de>'
__copyright__ = 'Copyright (C) 2020 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Prune angiosperm phylogenetic tree to family level and annotate ' \
		   'leaves with counts of avaible plastid genome entries in ' \
           'NCBI GenBank.'
__version__ = '2020.07.02.1700'


#############
# FUNCTIONS #
#############

def get_species_list_from_table(fp_table):
    '''
    Reads a csv file and returns the "organism" column as numpy array
    Params:
     - fp_table: file path to csv file
    '''
    plastome_table = pd.read_csv(fp_table, sep="\t", index_col=0)
    species_list = plastome_table["ORGANISM"].to_numpy()
    return species_list


def get_genome_count_per_family(species_list, species_in_family):
    '''
    Returns a dict with key = plant family name, value = number of sequenced plastid genomes for that family
    Params:
     - species_list: list of species names. One entry per sequenced plastid genome
     - species_in_family: dict with key=family name, value=set of species names
    '''
    genome_count_per_family = {}

    genome_count_per_species = Counter(species_list)

    for spec in genome_count_per_species:
        for family in species_in_family:
            if spec in species_in_family[family]:
                if family in genome_count_per_family:
                    genome_count_per_family[family] += genome_count_per_species[spec]
                else:
                    genome_count_per_family[family] = genome_count_per_species[spec]

    return genome_count_per_family


def get_family_from_lineage(lineage):
    '''
    Returns the family (taxonomic rank) name from a lineage
    Params:
     - lineage: dict with NCBI tax id as key and corresponding taxonomic name as value
    '''
    f_name = None
    for taxname in lineage.values():
        if taxname.endswith("ceae"):
            f_name = taxname
    return f_name


def get_species_in_family(species_set, ncbi_db):
    '''
    Constructs a dictionary from a set of species names where each family name key has a set of species name values assigned.
    Params:
     - species_set: the set of species to assign to families
     - ncbi_db: a NCBITaxa object to lookup a species' family name
    '''
    species_in_family = {}
    species_taxids = ncbi_db.get_name_translator(species_set)
    for species, taxid in species_taxids.items():
        f_name = get_family_from_lineage(ncbi_db.get_taxid_translator(ncbi_db.get_lineage(taxid[0])))
        if f_name:
            if f_name in species_in_family:
                species_in_family[f_name].add(species)
            else:
                species_in_family[f_name] = set([species])

    return species_in_family


def prune_to_family(tree, species_in_family):
    '''
    Prunes a tree from species-level leaves to family-level leaves
    Params:
     - tree: an ete3 Tree with species-level taxonomic names as leaf names
     - species_in_family: dict where each family name key has a set of species name values assigned
    '''

    # TODO: Check that detaching a child does not detach leaves of a different family
    # Probably an issue with some common ancestors being identified as close to the root
    #
    fam_nodes = []
    for family, species in species_in_family.items():
        if len(species) == 1:
            # if there is only on species in a family, simply rename the leaf
            fam_node = tree.get_leaves_by_name(list(species)[0].replace(" ", "_"))[0]
            fam_node.name = family
            fam_nodes.append(fam_node)
        else:
            fam_node = tree.get_common_ancestor([spec.replace(" ", "_") for spec in species])
            fam_node.name = family
            fam_nodes.append(fam_node)
    tree.prune([f_node.name for f_node in fam_nodes])


def get_species_count_per_family(species_in_family):
    '''
    Returns the number of species in each family
    Params:
     - species_in_family: dict where each family name key has a set of species name values assigned
    Returns:
     -  count_species_in_family: a dict where each family name key has the number of species in it assigned
    '''
    count_species_in_family = {}
    for k,v in species_in_family.items():
        count_species_in_family[k] = len(v)
    return count_species_in_family


def attach_counts_to_tree(tree, genome_count_per_family, species_count_per_family):
    '''
    Attaches TextFaces to tree leaves with plastid genome counts per leaf (=family) and leaf name
    Params:
     - tree: ete3 tree with family names as leaves
     - genome_count_per_family: dict with key = family name, value = number of sequenced plastid genomes for that family
    '''
    for leaf in tree.iter_leaves():
        leaf.add_face(TextFace(str(leaf.name)), column=0, position="aligned")
        if leaf.name in genome_count_per_family:
            leaf.add_face(TextFace("(" + str(species_count_per_family[leaf.name]) + "/" + str(genome_count_per_family[leaf.name]) + ")"), column=1, position="aligned")
        else:
            leaf.add_face(TextFace("(" + str(species_count_per_family[leaf.name]) + "/0)"), column=1, position="aligned")


def extend_leaf_branches(tree):
    '''
    -- DOES NOT YET WORK AS DESIRED --
    Extends leaf branches so the leaf nodes are aligned with the most distant leaf
    Params:
     - tree: ete3 tree
    '''
    maxDist = tree.get_farthest_leaf()[1]
    for leaf in tree.iter_leaves():
        dist = tree.get_distance(leaf)
        leaf.dist += maxDist - dist


def main(args):

    # STEP 1: Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)

    # STEP 2: Retrieve and/or update localized NCBI Taxonomy database
    ncbi = NCBITaxa()
    if (time.time() - os.path.getmtime(os.path.join(Path.home(), ".etetoolkit/taxa.sqlite"))) > 604800:
        ncbi.update_taxonomy_database()

    # STEP 3: Prune species-level tree to family-level

        # Step 3.1 Read tree from input file
    log.debug("Loading Tree...")
    t = Tree(args.infn, format=5)
        # STEP 3.2: Add species names to species_set_from_tree set
    log.debug("Gathering species(leaf) names...")
    species_set_from_tree = set()
    for leaf in t.iter_leaves():
        species_set_from_tree.add(leaf.name.replace("_"," "))
        # STEP 3.3: Assign species to families
    log.debug("Constructing dict of species in family...")
    species_in_family = get_species_in_family(species_set_from_tree, ncbi)
        # STEP 3.4: Prune the tree
    log.debug("Pruning Tree to family level...")
    prune_to_family(t, species_in_family)

    # STEP 4: Calculate counts of species per family and plastid genome entries per family and attach them to the tree leaves

        # STEP 4.1: Read plastid genome information from input table
    species_list_from_table = get_species_list_from_table(args.tablefn)
        # STEP 4.2: Count plastid genome entries per family
    log.debug("Counting plastid genome entries per family...")
    genome_count_per_family = get_genome_count_per_family(species_list_from_table, species_in_family)
        # STEP 4.3: Attach counts to tree leaves
    log.debug("Attaching counts to Tree...")
    attach_counts_to_tree(t, genome_count_per_family, get_species_count_per_family(species_in_family))

    # STEP 5: Set TreeStyle and render tree
    ts = TreeStyle()
    ts.mode = "c"
    ts.draw_guiding_lines = True
    ts.show_leaf_name = False
    log.debug("Rendering Tree...")
    t.render(args.outfn, w=10000, h=10000, tree_style=ts)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("-i", "--infn", type=str, required=True, help="path to input tree file. File has to be in Newick format 5 as specified in http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees")
    parser.add_argument("-o", "--outfn", type=str, required=True, help="path to .png output file")
    parser.add_argument("-t", "--tablefn", type=str, required=True, help="path to .csv file containing information on plastid genomes available in NCBI GenBank")
    args = parser.parse_args()
    main(args)
