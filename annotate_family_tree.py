from collections import Counter
import pandas as pd
import sys
import subprocess
import xml.etree.ElementTree as ET
from ete3 import Tree, TreeStyle, TextFace


class Lineage_Lookup:

    def __init__(self, fp_dump):
        print("Intializing lineage lookup...")
        self.lin_table = self.read_lineage_table(fp_dump)

    def read_lineage_table(self, fp_dump):
        '''
        Reads the fullnamelineage.dmp file that is part of the "new_taxdump" archive provided by NCBI at ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump
        Params:
         - fp_dump: file path to the .dmp file
        Returns:
         - lin_table: pandas DataFrame with columns SPECIES (index) and LINEAGE
        '''
        rem_tabs = lambda x : (x.replace("\t", "").strip())
        lin_table = pd.read_csv(fp_dump, sep="|", header=None, usecols=[1,2], names=["SPECIES", "LINEAGE"], converters = {"SPECIES": rem_tabs, "LINEAGE": rem_tabs}, index_col="SPECIES")
        return lin_table

    def reduce_to_spermatophyta(self):
        self.lin_table = self.lin_table.loc[self.lin_table["LINEAGE"].str.find("Spermatophyta") != -1]

    def get_family_name(self, species_name):
        '''
        Returns the taxonomic family name of a provided species name.
        Params:
         - species_name: str. scientific name of a species
        Returns:
         - family_name: str or None. family name of species_name
        '''
        family_name = None
        if species_name in self.lin_table.index:
            lin = self.lin_table.at[species_name, "LINEAGE"]
            for rank in lin.split(";"):
                if rank.endswith("ceae"):
                    family_name = rank.strip()
        return family_name


def get_species_list_from_table(fp_table):
    '''
    Reads a csv file and returns the "organism" column as numpy array
    Params:
     - fp_table: file path to csv file
    '''
    plastome_table = pd.read_csv(fp_table, sep="\t", index_col=0)
    species_list = plastome_table["ORGANISM"].to_numpy()
    return species_list

"""
def add_genome_counts(tree, species_list):
    '''
    Adds the number of sequenced plastid genomes per species to tree leaves.
    Params:
     - tree: The tree to edit
     - species_list: a pandas table containing plastome information
    '''
    species_list = plastome_table["ORGANISM"].to_numpy()
    counts_per_organism = Counter(species_list)
    for organism in counts_per_organism:
        for leaf in tree.get_leaves_by_name(organism.replace(" ", "_")):
            leaf.add_feature("sequenced_plastid_genomes", counts_per_organism[organism])
"""

def get_counts_per_family(species_list, species_in_family):
    '''
    Returns a dict with key = plant family name, value = number of sequenced plastid genomes for that family
    Params:
     - species_list: list of species names. One entry per sequenced plastid genome
     - species_in_family: dict with key=family name, value=set of species names
    '''
    counts_per_family = {}

    counts_per_species = Counter(species_list)

    for spec in counts_per_species:
        for family in species_in_family:
            for species in species_in_family[family]:
                if species == spec:
                    if family in counts_per_family:
                        counts_per_family[family] += counts_per_species[spec]
                    else:
                        counts_per_family[family] = counts_per_species[spec]

    return counts_per_family

"""
TM: turns out querying thousands of taxonomic entries takes a long while. Replaced this with a localized solution
def get_family_name(species_name):
    '''
    Returns the family name of a given species by querying GenBank.
    Params:
     - species_name: the scientific name of a plant species
    '''
    family_name = None

    species_name = species_name.replace("_", " ")

    esearchargs = ["esearch", "-db", "taxonomy", "-query", species_name]
    esearch = subprocess.Popen(esearchargs, stdout=subprocess.PIPE)
    efetchargs = ["efetch", "-mode", "xml"]
    efetch = subprocess.Popen(efetchargs, stdin=esearch.stdout, stdout=subprocess.PIPE)
    out, err = efetch.communicate()

    try:
        root = ET.fromstring(out)
        lineage = root.find("Taxon").find("LineageEx")
        for taxon in lineage.findall("Taxon"):
            if taxon.find("Rank").text == "family":
                family_name = taxon.find("ScientificName").text
                break
    except Exception as err:
        print("Couldn't find family name for species: " + str(species_name) + "\n" + str(err))

    return family_name
"""

def get_species_in_family(species_list, lineage_lookup):
    '''
    Constructs a dictionary from a list of species names where each family name key has a set of species name values assigned.
    Params:
     - species_list: the list of species to assign to families
    '''
    print("Executing get_species_in_family")
    species_in_family = {}
    for species in species_list:
        f_name = lineage_lookup.get_family_name(species)
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


def attach_counts_to_tree(tree, counts_per_family):
    '''
    Attaches TextFaces to tree leaves with plastid genome counts per leaf (=family) and leaf name
    Params:
     - tree: ete3 tree with family names as leaves
     - counts_per_family: dict with key = family name, value = number of sequenced plastid genomes for that family
    '''
    for leaf in tree.iter_leaves():
        leaf.add_face(TextFace(str(leaf.name)), column=0, position="aligned")
        if leaf.name in counts_per_family:
            leaf.add_face(TextFace(str(counts_per_family[leaf.name])), column=1, position="aligned")
        else:
            leaf.add_face(TextFace("0"), column=1, position="aligned")

# For some reason, the leaf nodes still do not align
def extend_leaf_branches(tree):
    '''
    Extends leaf branches so the leaf nodes are aligned with the most distant leaf
    Params:
     - tree: ete3 tree
    '''
    maxDist = tree.get_farthest_leaf()[1]
    for leaf in tree.iter_leaves():
        dist = tree.get_distance(leaf)
        leaf.dist += maxDist - dist

def color_node(node, color):

    for child in node.get_children():
        color_node(child)

def increment_color(color):


# TODO: Automatisiert die neueste fullnamelineage.dmp Datei vom NCBI ftp ziehen?
# -> ODER in ete3 eingebaute lokale NCBI Taxonomy DB nutzen
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Missing or too many arguments! Usage: some_functions.py [TREE_FILEPATH_IN] [TREE_FILEPATH_OUT] [TABLE_FILEPATH] [LINEAGE_DUMP]")
    else:
        species_list_table = get_species_list_from_table(sys.argv[3])
        species_list_tree = []
        print("Loading Tree...")
        t = Tree(sys.argv[1], format=5)
        print("Gathering species(leaf) names...")
        for leaf in t.iter_leaves():
            species_list_tree.append(leaf.name.replace("_"," "))
        print("Constructing lineage lookup...")
        lin_lookup = Lineage_Lookup(sys.argv[4])
        print("Reducing lineage lookup to Spermatophyta...")
        lin_lookup.reduce_to_spermatophyta()
        print("Constructing dict of species in family...")
        species_in_family = get_species_in_family(species_list_tree, lin_lookup)
        print("Counting species per family...")
        counts_per_family = get_counts_per_family(species_list_table, species_in_family)
        print("Pruning Tree to family level...")
        prune_to_family(t, species_in_family)
        print("Attaching counts to Tree...")
        attach_counts_to_tree(t, counts_per_family)
        print("Extending leaf branches...")
        extend_leaf_branches(t)

        ts = TreeStyle()
        ts.mode = "c"
        print("Rendering Tree...")
        t.render(sys.argv[2], w=10000, h=10000, tree_style=ts)
        #t.write(features=["sequenced_plastid_genomes"], outfile=sys.argv[2], format=5) # Keine Ahnung warum das Feature nicht mitgeschrieben wird
