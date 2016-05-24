#' fdrCsSAM
#' 
#' Estimates the false discovery rate for the identified cell-specific
#' differences in gene expression.
#' 
#' 
#' @param G Matrix of gene expression, columns ordered in the same order at the
#' cell-frequency matrix (n by p, n samples, p genes)
#' @param cc Matrix of cell-frequency. (n by k, n samples, k cell-types)
#' @param y A numeric vector of group association of each sample. Either 1 or 2.
#' @param n A nuermic vector describing the number of samples in a group
#' @param numcell The number of cell-types to consider
#' @param numgene The number of genes being considered
#' @param rhat The contrast in cell-type expression for each cell-type as
#' observed between the two groups being compared.
#' @param nperms The number of permutations to perform.
#' @param alternative Type of test to conduct - choose between
#' 'two.sided','greater',or 'less'
#' @param standardize Standardize sample or not. Default is TRUE
#' @param medianCenter Median center rhat distributions. Default is TRUE.
#' @param logRm Exponentiate data for deconvolution stage. Default is FALSE
#' @param logBase Base of logaritm used to determine exponentiation factor.
#' Default is 2
#' @param nonNeg For single channel arrays. Set any cell-specific expression
#' estimated as negative, to a ceiling of 0. It is conservative in its study of
#' differential expression. Default is FALSE.
#' @return A list.
#' \item{fdr.g}{A matirx false dicovery rates for csSAM comparison for each
#' cell-type at different thresholds. A set of 100 theresholds is determined
#' automatically from the data (k by 100, where k is number of cells).}
#' \item{avrhatperm}{}
#' \item{rhatperm}{A matrix sized pXkXg which stores the contrast of a
#' given gene g in cell type k in permutation p of the data.}
#' \item{cutp.g}{A matrix k by 100, where k is the number of cell tpes.
#' Lists the 100 cutoff thresholds for each cell-type as determined
#' automatically from the computed contrast.}
#' \item{rhat}{A matrix object with the result of contrasting the average
#' cell-specific expression profile of the two groups, per cell-type (Size k by
#' g where k is the number of cells and g is the number of genes).}
#' \item{ncall.g}{Number of genes called significant at the given cutoff
#' threshold with a FDR matching that indicated in fdr.g}
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @cite Shen-Orr2010
fdrCsSAM <-
function (G,cc,y,n,numcell,numgene,rhat,nperms,alternative='two.sided',standardize=TRUE,medianCenter=TRUE,logRm=FALSE,logBase = 2,nonNeg=FALSE) {
  numgene=ncol(G)
  rhatperm <- array(dim = c(nperms,numcell,numgene))
  perm = list()
  
  for (i in 1:nperms) {
    o=sample(1:length(y))
    ystar = y[o]
    for (curset in 1:2) {
      perm[[curset]]= csfit(cc[ystar==curset,], G[ystar==curset,],logRm,logBase)
    }
    rhatperm[i,,] = csSAM(perm[[1]]$ghat, perm[[1]]$se, n[1], perm[[2]]$ghat, perm[[2]]$se, n[2],
              standardize, medianCenter, nonNeg)
  }
  cutp.g=matrix(NA,nrow=numcell,ncol=100)
  numcut = ncol(cutp.g)
  
  fdr.g=ncall.g=nperm.g<-array(dim = c(numcell, numcut))
  
  for(j in 1:numcell)
    cutp.g[j,]=seq(0,max(abs(rhat[j,])),length=100)	
  for (i in 1:numcut) {
    for (curcell in 1:numcell) {
		if(alternative == 'two.sided') {
			fdr.g[curcell,i]=sum(abs(rhatperm[,curcell,])>cutp.g[curcell,i])/nperms /sum(abs(rhat[curcell,])>cutp.g[curcell,i])
			ncall.g[curcell,i]=sum(abs(rhat[curcell,])>cutp.g[curcell,i])
		}
		if(alternative == 'greater') {
			fdr.g[curcell,i]=sum(rhatperm[,curcell,]>cutp.g[curcell,i])/nperms /sum(rhat[curcell,]>cutp.g[curcell,i])
			ncall.g[curcell,i]=sum(rhat[curcell,]>cutp.g[curcell,i])
		}
		if(alternative == 'less') {
			# [RG] BUG FIX: should be < - cutp.g[curcell,i]
#			fdr.g[curcell,i]=sum(rhatperm[,curcell,]< -cutp.g[curcell,i])/nperms /sum(rhat[curcell,]>cutp.g[curcell,i])
			fdr.g[curcell,i]=sum(rhatperm[,curcell,]< -cutp.g[curcell,i])/nperms /sum(rhat[curcell,] < -cutp.g[curcell,i])
			ncall.g[curcell,i]=sum(rhat[curcell,]< -cutp.g[curcell,i])
		}		
    }
  }
  
  fdr.g =pmin(fdr.g,1)
  
  for (j in 1:numcell)	 {
    fdr.g[j,]=make.monotone(fdr.g[j,])
  }	
  
  return (list(fdr.g=fdr.g, rhatperm = rhatperm, cutp.g = cutp.g, rhat = rhat,ncall.g = ncall.g, alternative))
}

