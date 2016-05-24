#' findSigGene
#' 
#' Find the false discovery rate for each gene in each cell-type.
#' 
#' 
#' @param G Gene expression matrix of heterogenous tissue measurements
#' @param cc Matrix of cell-frequency measures per person
#' @param y Numeric group association of each sample. Either 1 or 2.
#' @param rhat Matrix of cell-specific contrasts for each gene in each cell-type
#' as computed for the original group classification.
#' @param csSAMData List object returned from fdrCsSAM.
#' @return A matrix size k by g where k is the number of cell-types and g is the
#' number of genes. For each cell in the matirx, listed is the FDR of the gene
#' for a difference in a given cell-type.
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @cite Shen-Orr2010
findSigGene <- function(G,cc,y,rhat,csSAMData) {
# @useDynLib csSAM
#	.Call('findSigGenes', rhat, csSAMData$cutp.g, csSAMData$fdr.g, PACKAGE='csSAM')
# old plain R version
#.findSigGene <- function(G,cc,y,rhat,csSAMData) {
  numgene=ncol(G)
  numcell = ncol(cc)
  thresholdVec = csSAMData$fdr.g
  cutoff <- csSAMData$cutp.g
  thresholdLen = length(thresholdVec[numcell,])
  sigGene <- array(dim = c(numcell, numgene))
  sigGene[,] = 1
  
  for (curThresh in 1:thresholdLen) {
    for (curcell in 1:numcell) {
      for (curgene in 1:numgene) {
        if(abs(rhat[curcell,curgene]) >= abs(cutoff[curcell,curThresh])) {
          sigGene[curcell,curgene] = thresholdVec[curcell,curThresh]
        }
      }
    }
  }
  
  return (sigGene)
}

findSigGene <- compiler::cmpfun(findSigGene) 
