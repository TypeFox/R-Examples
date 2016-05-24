##########################################################
# Recode a genotype of n genes with k<=n *distinct*      #
# alleles 1, 2, ..., k in such a way that the *distinct* #
# alleles of the recoded genotype appear in the natural  #
# order 1, 2, ..., k.                                    #
# For shortness, I'll say that any recoded genotype is   #
# "ordered".                                             #
# NOTE: If perm is already ordered, do nothing.          #
# Author: Diego Zardetto
##########################################################
`recode` <- function(perm){
  gr <- unique(perm)
  ngr <- length(gr)
  gr.rec <- 1:ngr
  if (identical(gr, gr.rec)) return(perm)
  perm.rec <- factor(perm, levels=gr)
  levels(perm.rec) <- gr.rec
  as.integer(unclass(perm.rec))
}