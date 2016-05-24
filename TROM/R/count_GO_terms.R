count_GO_terms <-
function(gene_lists, gene_pop, GO_pop, GO2geneID) {
  sapply(gene_lists, FUN=function(genes) {
    n <- length(genes) # num of genes in a sample
    N <- length(gene_pop) # num of genes in the population
    sapply(GO_pop, FUN=function(GO_ID) {
      spec_genes <- GO2geneID[[GO_ID]]
      M <- sum(spec_genes %in% gene_pop) # num of genes that have GO terms
      m_obs <- sum(spec_genes %in% genes) # num of stage-associated genes that have GO terms
      p_val <- phyper(q=m_obs, m=M, n=N-M, k=n, lower.tail=FALSE) + dhyper(x=m_obs, m=M, n=N-M, k=n)
      return(p_val)    	
    })
  })
}
