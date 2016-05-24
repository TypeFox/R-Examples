# Changes:
#	2013-09-19: Adjust the input to the new getSigTest output
#	2013-09-21: Still some adjustments were necessary...

plot.re <- function(x,...){
  
  if(x$inputN>1){
     pvalues <- x[[1]]$pvalues
     rejLine <- x[[1]]$method
     alpha <- x[[1]]$alpha
     
      for(i in 2:x$inputN){
        pvalues <- rbind(pvalues,x[[i]]$pvalues)
      }
  } else {
      pvalues <- x$pvalues
      rejLine <- x$method
      alpha <- x$alpha
  }

  possibleRL <- c("bh","bonferroni","simes")
  if(sum(is.element(rejLine,possibleRL))==0){
    rejLine <- NULL
    alpha <- NULL
  }
  
  rejectionPlot(X=pvalues, rejLine=rejLine, alpha=alpha, ...)
}