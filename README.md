# KEGG_edge_list
Make edge list from KEGG KGML files

This function will compute an undirected interaction edge list from the KEGG pathway files which can be used as input for the CPI interaction prediction tool BiPredict.

## Requirements
FIRST: create a new directory and move the file 'kegg_pathway_network.R' and 'kegg_kgml.sh' to this directory.
You need this script and the file 'kegg_kgml.sh' in the same directory. Make sure the file 'kegg_kgml.sh' is executable (chmod u+x kegg_kgml.sh).
During the process, the KEGG KGML files will be downloaded and stored in a directory which also will be created and named 'KGML_files'. If you wish to keep the KGML files, chose 'delete_dir = FALSE', otherwise they will be deleted after processing.
MANDATORY: you need to delete or move the directory 'KGML_files' after each species processing (easy way: use 'delete_dir = TRUE')!!!


## Parameters
```
ORGANISM:  KEGG organism abbreviation, please find them all here: https://www.genome.jp/kegg/catalog/org_list.html
delete_dir: If the folder with downloaded KEGG KGML files should be deleted after processing
only_c_compounds: If only compounds with a compound 'C' number in KEGG should be included in the edge list (recommended)
```
## Usage 
```
start R
source 'source kegg_pathway_network.R'
edge_list = make_edge_list_from_kegg(ORGANISM,delete_dir=T,only_c_compounds=T)
edge_list_ath = make_edge_list_from_kegg('ath',delete_dir=T,only_c_compounds=T)
```
The resulting edge list will be saved to a file named 'kegg_edge_list_ORGANISM.csv'
