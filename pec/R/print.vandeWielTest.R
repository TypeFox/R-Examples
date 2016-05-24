##' @export
print.vandeWielTest <- function(x,eps=0.0001,pdigits=4,...){
  cat("\nvan de Wiel test based on ",x$B," data splits\n")
  cat("\nTraining sample size: ",x$M,"\n")
  cat("\nTest sample size: ",x$N-x$M,"\n")
  if (length(x$testIBS)==2){
    cat("\nP-values based on integrated Brier score residuals:")
    cat("\nRange of integration: [",x$testIBS[1],"--",x$testIBS[2],"]\n\n")
    ibsP <- sapply(x$Comparisons,function(x)x$pValueIBS)
    ibsP <- format.pval(ibsP,digits=pdigits,eps=eps)
    ibsMat <- matrix(ibsP,ncol=1)
    rownames(ibsMat) <- names(x$Comparisons)
    colnames(ibsMat) <- "p-value (IBS)"
    print(ibsMat,quote=FALSE,...)
  }
  NT <- length(x$testTimes)
  if (NT>0){
    cat("\nMatrix of time point wise p-values:\n\n")
    if (NT>5){
      showTimes <- sort(sample(x$testTimes))
      showTimePos <- prodlim::sindex(jump.times=x$testTimes,eval.times=showTimes)
    }
    else{
      showTimes <- x$testTimes
      showTimePos <- 1:NT
    }
    mat <- do.call("rbind",lapply(x$Comparisons,function(comp){
      format.pval(comp$pValueTimes[showTimePos],digits=pdigits,eps=eps)
    }))
    colnames(mat) <- paste("t=",showTimes)
    print(mat,quote=FALSE,...)
  }
  invisible(mat)
}
