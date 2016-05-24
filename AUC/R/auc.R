#' Compute the area under the curve of a given performance measure.
#'
#' This function computes the area under the sensitivity curve (AUSEC), the area under the specificity curve (AUSPC), 
#' the area under the accuracy curve (AUACC), or the area under the receiver operating characteristic curve (AUROC).
#' 
#' @param x an object produced by one of the functions \code{sensitivity}, \code{specificity},  \code{accuracy}, or \code{roc}
#' @param min a numeric value between 0 and 1, denoting the cutoff that defines the start of the area under the curve
#' @param max a numeric value between 0 and 1, denoting the cutoff that defines the end of the area under the curve
#' 
#' @examples
#' 
#' data(churn)
#' 
#' auc(sensitivity(churn$predictions,churn$labels))
#' 
#' auc(specificity(churn$predictions,churn$labels))
#' 
#' auc(accuracy(churn$predictions,churn$labels))
#' 
#' auc(roc(churn$predictions,churn$labels))
#' 
#' 
#' @references Ballings, M., Van den Poel, D., Threshold Independent Performance Measures for Probabilistic Classifcation Algorithms, Forthcoming.
#' @seealso \code{\link{sensitivity}}, \code{\link{specificity}}, \code{\link{accuracy}}, \code{\link{roc}}, \code{\link{auc}}, \code{\link{plot}}
#' @return A numeric value between zero and one denoting the area under the curve
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@UGent.be}
auc <-  function(x, min=0,max=1) {
  
  if (any(class(x) == "roc")) {
        if (min != 0 || max != 1 ) {
          x$fpr <- x$fpr[x$cutoffs >= min & x$cutoffs <= max]
          x$tpr <- x$tpr[x$cutoffs >= min & x$cutoffs <= max]
        }
        
        
        ans <- 0
        for (i in 2:length(x$fpr)) {
          ans <- ans + 0.5 * abs(x$fpr[i] - x$fpr[i-1]) * (x$tpr[i] + x$tpr[i-1])
        }
        
  }else if (any(class(x) %in% c("accuracy","sensitivity","specificity"))) {
  
        if (min != 0 || max != 1 ) {
          x$cutoffs <- x$cutoffs[x$cutoffs >= min & x$cutoffs <= max]
          x$measure <- x$measure[x$cutoffs >= min & x$cutoffs <= max]
        }
        
        ans <- 0
        for (i in 2:(length(x$cutoffs))) {
          ans <- ans + 0.5 * abs(x$cutoffs[i-1] - x$cutoffs[i]) * (x$measure[i] + x$measure[i-1])
        }
            
  }
  
  
  
  return(as.numeric(ans))
}