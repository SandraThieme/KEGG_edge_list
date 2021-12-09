# KEGG_edge_list
Make edge list from KEGG KGML files

This function will compute an undirected interaction edge list from the KEGG pathway files which can be used as input for the CPI interaction prediction tool BiPredict.
#' FIRST: create a new directory and move the file 'kegg_pathway_network.R' and 'kegg_kgml.sh' to this directory
#' you need this script and the file 'kegg_kgml.sh' in the same directory
#' make sure the file 'kegg_kgml.sh' is executable (chmod u+x kegg_kgml.sh)
#' during the process, the KEGG KGML files will be downloaded and stored in a directory which also will be created and named 'KGML_files'
#' if you wish to keep the KGML files, chose 'delete_dir = FALSE', otherwise they will be deleted after processing
#' MANDATORY: you need to delete or move the directory 'KGML_files' after each species processing (easy way: use 'delete_dir = TRUE')!!!
#' start R
#' source the script: 'source kegg_pathway_network.R'
#' the resulting edge list will be saved to a file named 'kegg_edge_list_ORGANISM.csv'
#'
#' @param ORGANISM KEGG organism abbreviation, please find them all here: https://www.genome.jp/kegg/catalog/org_list.html
#' @param delete_dir If the folder with downloaded KEGG KGML files should be deleted after processing
#' @param only_c_compounds If only compounds with a compound 'C' number in KEGG should be included in the edge list (recommended)
#'
#' @example edge_list = make_edge_list_from_kegg(ORGANISM,delete_dir=T,only_c_compounds=T)
#' @example edge_list_ath = make_edge_list_from_kegg('ath',delete_dir=T,only_c_compounds=T)
