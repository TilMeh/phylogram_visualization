# phylogram_visualization

Create a family-level phylogram based on the Dated angiosperm phylogram file provided in (https://bdj.pensoft.net/article/39677/instance/5312825/ Suppl. material 5) using ETE Toolkit

Manual preparations:
The source phylogram is provided in a format that is not readable by ETE Toolkit. The following steps were needed to convert it to a readable Newick format:
* Extract the actual tree:
  $ awk '/^tree /' oo_330891.tre > largetree.tre
* Manually remove the leading and closing elements that do not belong to Newick format
* remove additional node attributes:
  $ sed "s/[[][^]]*[]]//g" largetree.tre > largetree_pruned.tre
* ... (explain further steps)
