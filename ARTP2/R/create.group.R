
create.group <- function(pathway, rs){
  
  pd <- pathway[pathway$SNP %in% rs, ]
  GeneInGroup <- unique(pd$Gene)
  ngene <- length(GeneInGroup)
  GeneIdx <- list()
  vGeneIdx <- NULL
  GeneStartEnd <- matrix(NA, ngene, 2)
  rownames(GeneStartEnd) <- GeneInGroup
  colnames(GeneStartEnd) <- c("Start", "End")
  N.SNP <- NULL
  for(g in 1:ngene){
    gene <- GeneInGroup[g]
    snps <- pd$SNP[pd$Gene == gene]
    N.SNP <- c(N.SNP, length(snps))
    GeneIdx[[g]] <- rep(NA, length(snps))
    for(i in 1:length(snps)){
      GeneIdx[[g]][i] <- which(rs == snps[i])
    }
    GeneStartEnd[gene, "Start"] <- length(vGeneIdx) + 1
    vGeneIdx <- c(vGeneIdx, GeneIdx[[g]])
    GeneStartEnd[gene, "End"] <- length(vGeneIdx)
  }
  
  list(GeneInGroup = GeneInGroup, GeneIdx = GeneIdx, 
       vGeneIdx = vGeneIdx, GeneStartEnd = GeneStartEnd, 
       N.SNP = N.SNP)
  
}


