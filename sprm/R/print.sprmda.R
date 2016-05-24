print.sprmda <-
function(x,...){
#  cat(paste(attr(x$coefficients, "Call"), "\n", sep=""))
  cat("Sparse partial M-robust regression \n")
  cat(paste(" Number of components: ", x$inputs$a))
  cat(paste("\n Sparsity parameter: eta=", x$inputs$eta))
  cat(paste(c(paste("\n weight function: ", x$inputs$fun, "with cutoff"), x$inputs$constants), collapse=" "))
}
