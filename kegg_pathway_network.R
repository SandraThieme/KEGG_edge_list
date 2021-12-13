## ----------------------------------------------------------------
##
##   R package
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
##   02110-1301 USA
##
## -----------------------------------------------------------------

#' @title Organism specific CPI edge list based on the KEGG pathway network
#' @author Sandra Thieme
#'
#' @description
#' This function will compute an undirected interaction edge list from the KEGG pathway files which can be used as input for the CPI interaction prediction tool BiPredict.
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



#get files from KEGG DB
download_kegg_kgml = function(reference){
  download.file(paste('http://rest.kegg.jp/list/pathway/',reference,sep = ''),paste('pathway_info_',reference,'.txt',sep = ''),method='wget')
  all_pathways = read.table(paste('pathway_info_',reference,'.txt',sep = ''),sep = '\t')
  all_pathways$V1 = gsub('path:','',all_pathways$V1)
  all_pathways$V3 = all_pathways$V1
  all_pathways$V3 = gsub('[a-z]','',all_pathways$V3)
  colnames(all_pathways) = c('pathway','name','no')
  hsa_pw_numbers = unique(all_pathways$no)
  for (pw in hsa_pw_numbers){
    print(pw)
    download.file(paste('http://rest.kegg.jp/get/',reference,pw,'/kgml',sep = ''),paste(reference,pw,'.xml',sep = ''),method='wget')
  }
}

make_edge_list_entry = function(PATH='./',NAME='entry'){
  reaction_files = list.files(PATH,pattern = NAME)
  edge_list = data.frame()
  for (k in 1:length(reaction_files)){
    print(reaction_files[k])
    ### read files with reactions after pre processing by kegg_kgml.sh ###
    entry = read.table(paste(PATH,reaction_files[k],sep = ''),sep = ';',fill = T,col.names = c('id','name','type','reaction'))
    edge_list= rbind(edge_list,entry)
  }
  return(edge_list)
}

make_edge_list_reaction = function(PATH='./',NAME='reaction'){
  reaction_files = list.files(PATH,pattern = NAME)
  edge_list = data.frame()
  for (k in 1:length(reaction_files)){
    print(reaction_files[k])
    ### read files with reactions after pre processing by kegg_kgml.sh ###
    reaction = read.table(paste(PATH,reaction_files[k],sep = ''),sep = ';',fill = T, col.names = c('V1','id','name','type','substrate.id','substrate.name','product.id','product.name'))
    reaction$V1<-NULL
    edge_list= rbind(edge_list,reaction)
  }
  return(edge_list)
}

#make one entry for each reaction
seperate_reactions = function(reaction){
  reaction$name = as.character(reaction$name)
  for (i in 1:dim(reaction)[1]){
    reaction_id = unlist(strsplit(reaction[i,'name'],';',fixed = T))
    if (length(reaction_id)>1){
      for (j in 1:(length(reaction_id)-1)){
        reaction[i,'name'] = reaction_id[1]
        last=dim(reaction)[1]+1
        reaction[last,]=reaction[i,]
        reaction[last,'name']=reaction_id[j+1]
      }
    }
  }
  return(reaction)
}

#make one entry for each reaction
seperate_entrys = function(entry){
  entry$reaction = as.character(entry$reaction)
  for (i in 1:dim(entry)[1]){
    entry_id = unlist(strsplit(entry[i,'reaction'],';',fixed = T))
    if (length(entry_id)>1){
      for (j in 1:(length(entry_id)-1)){
        entry[i,'reaction'] = entry_id[1]
        last=dim(entry)[1]+1
        entry[last,]=entry[i,]
        entry[last,'reaction']=entry_id[j+1]
      }
    }
  }
  entry$name = as.character(entry$name)
  for (i in 1:dim(entry)[1]){
    entry_id = unlist(strsplit(entry[i,'name'],';',fixed = T))
    if (length(entry_id)>1){
      for (j in 1:(length(entry_id)-1)){
        entry[i,'name'] = entry_id[1]
        last=dim(entry)[1]+1
        entry[last,]=entry[i,]
        entry[last,'name']=entry_id[j+1]
      }
    }
  }
  return(entry)
}

make_edge_list_from_kegg = function(ORGANISM,delete_dir=T,only_c_compounds=T){
  #make working directory
  system2('mkdir',args = ('KGML_files'))
  #set as working directory
  setwd("./KGML_files")
  #WD = getwd()

  download_kegg_kgml(ORGANISM)
  #do pre-processing using kegg_kgml.sh in bash
  system2('cp',args = c('../kegg_kgml.sh','.'))
  system2('./kegg_kgml.sh')
  entry_df = make_edge_list_entry()
  entry_df = entry_df[!is.na(entry_df$reaction),]
  reaction_df = make_edge_list_reaction()
  reaction_df = seperate_reactions(reaction_df)
  entry_df = seperate_entrys(entry_df)
  entry_df = entry_df[entry_df$type == 'gene',]
  rownames(entry_df)<-NULL
  entry_df$type <- NULL
  colnames(reaction_df)[2] = 'reaction'
  all_reactions_kegg = merge(reaction_df,entry_df)
  all_reactions_kegg = unique(all_reactions_kegg)
  all_reactions_kegg$substrate.name = gsub('cpd:','',all_reactions_kegg$substrate.name)
  all_reactions_kegg$product.name = gsub('cpd:','',all_reactions_kegg$product.name)
  sub_enz_hsa = unique(all_reactions_kegg[,c("substrate.name","name")])
  prod_enz_hsa = unique(all_reactions_kegg[,c("product.name","name")])
  colnames(prod_enz_hsa) = c('chemical','ProteinID')
  colnames(sub_enz_hsa) = c('chemical','ProteinID')
  kegg_edge_list = unique(rbind(sub_enz_hsa,prod_enz_hsa))
  kegg_edge_list$ProteinID = gsub(paste(ORGANISM,':',sep = ''),'',kegg_edge_list$ProteinID)
  rownames(kegg_edge_list)<-NULL

  if (only_c_compounds){
    print(c('length of edge list: ',nrow(kegg_edge_list)))
    kegg_edge_list = kegg_edge_list[grep('C',kegg_edge_list$chemical),]

    list1 = lapply(kegg_edge_list$chemical, FUN = function(x) unlist(strsplit(x,';',fixed = T)))
    list2 = lapply(list1, FUN = function(x) x[grep('C',x)])
    list2 = lapply(list2, FUN = function(x) paste(x,collapse = ';'))
    kegg_edge_list$chemical = unlist(list2)

    print(c('length of edge list with only compounds containing C number : ',nrow(kegg_edge_list)))
    print('If there is a huge difference, try: only_c_compounds=FALSE and take a look on the chemical IDs')
  }

  write.csv(kegg_edge_list,file = paste('../kegg_edge_list_',ORGANISM,'.csv',sep = ''),row.names = F)
  setwd("../")

  if(delete_dir){system2('rm',args = c('-r','KGML_files'))}

  rownames(kegg_edge_list)<-NULL
  return(kegg_edge_list)
}
