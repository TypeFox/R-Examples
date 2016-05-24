`summary.deseasonalize` <-
function(object, ...){
a <- object$dspar
if (nrow(a)==1) out <- a else {
  PL <- round(100*exp(-0.5*(a[,4]-min(a[,4]))),1)
  i <- order(PL, decreasing=TRUE)
  ind <- i[PL[i]>1]
  out <- cbind(a[ind,,drop=FALSE], PL[ind])
  colnames(out) <- c(colnames(a), "Plausibility %")
}
out
}

