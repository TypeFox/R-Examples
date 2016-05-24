write_top_GO_terms <-
function(output_file, associated_gene_list, all_genes, all_genes_w_GO, geneID2GO, topNum) {
  write("", file=output_file)
  top_GO_all <- sapply(1:length(associated_gene_list), FUN=function(i) {
    temp <- find_top_GO_terms(associated_gene_list[[i]], all_genes, all_genes_w_GO, geneID2GO, topNum)
    write(paste("Sample:", names(associated_gene_list)[i]), file=output_file, append=TRUE)
    write.table(temp, file=output_file, append=TRUE, quote=FALSE, sep="\t")
    write("\n", file=output_file, append=TRUE)
    return(temp)
  })
  return(top_GO_all)
}
