`print.homals` <-
function(x, ...)
{
  nvar <- dim(x$dframe)[2]
  cat("\nCall: ")
  print(x$call)
  cat("\nLoss:",x$loss,"\n")
  cat("\nEigenvalues:\n")
  eigen.val <- round(x$eigenvalues, 4)
  names(eigen.val) <- paste("D",1:x$ndim,sep="")
  print(eigen.val)
  cat("\n")
  
  loadmat <- t(sapply(x$loadings, function(xx) (xx[1,])))
  if (x$ndim == 1) loadmat <- t(loadmat)
  colnames(loadmat) <- paste("D", 1:x$ndim, sep = "")
  cat("Variable Loadings:\n")
  print(loadmat)
  cat("\n")
  invisible(x)
}

