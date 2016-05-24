summary.sprm <-
function(object,...){
  cat("Sparse partial M-robust regression \n")
  cat(paste(" Number of components: ", object$inputs$a))
  cat(paste("\n Sparsity parameter: eta =", object$inputs$eta))
  cat(paste(c(paste("\n weight function: ", object$inputs$fun, "with cutoff"), object$inputs$constants), collapse=" "))
  cat("\n\nNumber of variables included in the model: ")
  #  print(paste(head(object$coefficients)))
  cat(paste0("\n With ", 1, " component : "))
  cat(length(object$used.var[[2]]))
  if(object$inputs$a>1){
    for(i in 2:object$inputs$a){
      cat(paste0("\n With ", i, " components: "))
      cat(length(object$used.var[[2*i]]))
    }
  }
  cat("\n\nPercentage of explained variance\n")
  Vartab <- data.frame(X=cumsum(object$Xvar), y=cumsum(object$Yvar))
  rownames(Vartab) <- c(paste("With ", 1:object$inputs$a, " component(s):"))
  print(Vartab,...)
  
#  cat("\n\nResiduales\n")
#  cat(summary(object$residuals))
}
