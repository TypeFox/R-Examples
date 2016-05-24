#' runSAM
#' 
#' A lightweight version of the SAM algorithm, only performs two group
#' comparison with equal deltas on each tail
#' 
#' 
#' @param G Matrix of gene expression, columns ordered in the same order at the
#' cell-frequency matrix (n by p, n samples, p genes)
#' @param y Numeric group association of each sample. Either 1 or 2.
#' @param s0.sam Input or computed value of SAM exchangeability factor. Default
#' is determined automatically
#' @param stand.r Median center and standardize arrays. Default is TRUE.
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @cite Shen-Orr2010
runSAM <-
function(G,y,s0.sam=NULL,stand.r=TRUE) {
  
  if(stand.r == TRUE) {
	G=scale(t(G),center=apply(G,1,median),scale=F)
	}
  if(is.null(s0.sam)){
	s0.sam=quantile(ttest.func(G,y)$sd,.5,na.rm=TRUE)
	}
  tt.sam=ttest.func(G,y,s0=s0.sam)$tt
  return (tt.sam)
}

