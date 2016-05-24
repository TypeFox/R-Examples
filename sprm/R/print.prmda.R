print.prmda <-
function(x,...){
  cat("Partial M-robust regression \n")
  cat(paste(" Number of components: ", x$inputs$a))
  cat(paste(c(paste("\n weight function: ", x$inputs$fun, "with cutoff"), x$inputs$constants), collapse=" "))
}
