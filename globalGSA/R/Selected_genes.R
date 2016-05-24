Selected_genes <-
function(Gene, gene_list, data) {
  idX <- as.character(gene_list$Id[gene_list$Gene == Gene])
  return(data[,idX])
}
