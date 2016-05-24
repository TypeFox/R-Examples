summary.prmda <-
function(object,...){
  cat("Partial M-robust regression \n")
  cat(paste(" Number of components: ", object$inputs$a))
  cat(paste(c(paste("\n weight function: ", object$inputs$fun, "with cutoff"), object$inputs$constants), collapse=" "))
  cat("\n\nPercentage of explained variance\n")
  Vartab <- data.frame(X=cumsum(object$Xvar), y=cumsum(object$Yvar))
  rownames(Vartab) <- c(paste("With ", 1:object$inputs$a, " component(s):"))
  print(Vartab, ...)
  
}
