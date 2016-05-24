#' csSAM
#' 
#' Computes the constrast between groups for the deconvolved cell-specific
#' expression for each cell-type
#' 
#' 
#' @param ghat1 Expression matrix of deconvolved cell-specific gene expression
#' estimates for group 1.
#' @param se1 Standard error group 1
#' @param n1 Group 1 size
#' @param ghat2 Expression matrix of deconvolved cell-specific gene expression
#' estimates for group 2.
#' @param se2 Standard error group 2
#' @param n2 Group 2 size
#' @param standardize Standardize contrast values
#' @param medianCenter Median center rhat distributions for each cell-type
#' @param nonNeg Negative values not allowed such as in a single channel
#' microarray. Zero them if negative (a conervative option)
#' @return A matrix object with the result of contrasting the average
#' cell-specific expression profile of the two groups, per cell-type (Size k by
#' g where k is the number of cells and g is the number of genes).
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @cite Shen-Orr2010
csSAM <-
function(ghat1, se1, n1, ghat2, se2, n2, standardize,
                  medianCenter=TRUE, nonNeg=FALSE) {
  numcell=nrow(ghat1)
  if (nonNeg) {
    ghat1[ghat1 < 0] = 0	
    ghat2[ghat2 < 0] = 0
  }
  rhat = ghat2 - ghat1
  
  ##if the data is to be standardized
  if (standardize) {
    se=((n1*se1^2+n2*se2^2)/(n1+n2))^(1/2)
    s0.r=apply(se,1,quantile, .5,na.rm=TRUE)
    for (i in 1:numcell) {
      rhat[i,]=rhat[i,]/(se[i,]+s0.r[i])
    }
  }
  
  ##if the arrays are to be median-centered
  if (medianCenter) {
    for(i in 1:numcell) {
      rhat[i,]=rhat[i,]-median(rhat[i,])
    }
  }
  return (rhat)
}

