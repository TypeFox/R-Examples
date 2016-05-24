#' fdrSAM
#' 
#' Calculate the false discovery rate (FDR) by permutation for the group
#' differences as calculated by SAM.
#' 
#' 
#' @param G Matrix of gene expression, columns ordered in the same order at the
#' cell-frequency matrix (n by g, n samples, g genes)
#' @param y A numeric vector of group association of each sample. Either 1 or 2.
#' @param nperms Number of permutations to run. User responsability to the
#' number appropriately fitting the sample size.
#' @param tt.sam Real group comparison t-test statistic value
#' @param alternative Type of test. Choices are 'two.sided','greater' or 'less'
#' @return A list
#' \item{fdr.sam}{A vector false dicovery rates for SAM comparison at
#' different thresholds. A set of 100 theresholds is determined automatically
#' from the data.}
#' \item{ncall.sam}{Number of genes called significant at the given cutoff
#' threshold with a FDR matching that indicated in fdr.sam}
#' \item{ttstar.sam}{A matrix listing the t statistic for each gene in each
#' permutation. (p by g, p permutations, g genes)}
#' \item{sigGene.sam}{A vector of length equal to the number of genes being
#' considered. For each gene the estimated FDR is listed.}
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @cite Shen-Orr2010
fdrSAM <-function (G,y,nperms,tt.sam,alternative= 'two.sided') {
  numgene=ncol(G)
  perm = list()
  ttstar.sam <- array(dim = c(nperms,numgene))

  
  for (i in 1:nperms) {
    o=sample(1:length(y))
    ystar = y[o]
    ttstar.sam[i,] = runSAM(G, ystar)
  }
  cutp.sam=seq(0,max(abs(tt.sam)),length=100)
  sigGene.sam <- vector(mode = 'numeric',numgene)
  
  fdr.sam = ncall.sam = rep(NA, 100)
  

  for(i in 1:100) {
	if(alternative == 'two.sided') {
		fdr.sam[i]=(sum(abs(ttstar.sam)>cutp.sam[i])/nperms)/sum(abs(tt.sam)>cutp.sam[i])
		ncall.sam[i]=sum(abs(tt.sam)>cutp.sam[i])
		sigGene.sam[which(abs(tt.sam)>cutp.sam[i])] = fdr.sam[i]
	}
	if(alternative == 'greater') {
		fdr.sam[i]=(sum(ttstar.sam>cutp.sam[i])/nperms)/sum(tt.sam>cutp.sam[i])        
		ncall.sam[i]=sum(tt.sam>cutp.sam[i])
		sigGene.sam[which(tt.sam>cutp.sam[i])] = fdr.sam[i] 
	}
	if(alternative == 'less') {
      fdr.sam[i]=(sum(ttstar.sam< -cutp.sam[i])/nperms)/sum(tt.sam< -cutp.sam[i])   
      ncall.sam[i]=sum(tt.sam < -cutp.sam[i])
	  sigGene.sam[which(tt.sam< -cutp.sam[i])] = fdr.sam[i] 
    }
  }
  
  fdr.sam=pmin(fdr.sam,1)
  fdr.sam=make.monotone(fdr.sam)
  
  return (list(fdr.sam = fdr.sam, ncall.sam = ncall.sam,ttstar.sam = ttstar.sam,sigGene.sam = sigGene.sam))
}

