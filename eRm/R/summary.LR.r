summary.LR <- function(object,...)
# summary method for objects of class "LR" (from LRtest")
{
  cat("\n")
  cat("Andersen LR-test: \n")
  cat("LR-value:", round(object$LR,3),"\n")
  cat("Chi-square df:",object$df,"\n")
  cat("p-value: ",round(object$pvalue,3),"\n")
  cat("\n")

  mt_vek <- apply(object$X,2,max,na.rm=TRUE)

  for (i in 1:length(object$betalist)) {
    cat("\n")
    cat("Subject Subgroup: ",object$spl.gr[i],":",sep="")
    cat("\n")
    cat("Log-likelihood: ",object$likgroup[i])
    cat("\n\n")
    cat("Beta Parameters: \n")
    betavec <- object$betalist[[i]]
    if (!all(is.na(object$selist[[i]]))) {
      coeftable <- rbind(betavec,object$selist[[i]])
      rownames(coeftable) <- c("Estimate","Std.Err.")
      print(coeftable)
    } else {
      print(betavec)
    }
    cat("\n")
  }
}


