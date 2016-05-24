find_top_GO_terms <-
function(associated_genes, all_genes, all_genes_w_GO, geneID2GO, topNum) {
  myInterestingGenes <- intersect(associated_genes, all_genes_w_GO)
  geneNames <- intersect(all_genes, all_genes_w_GO)
  enrichment_test(myInterestingGenes, geneNames, geneID2GO, topNum)
}
