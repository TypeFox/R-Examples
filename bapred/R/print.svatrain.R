print.svatrain <-
function(x, ...) {

  cat("'SVA'-adjusted training data with information for addon batch effect adjustment.", "\n")
  cat(paste("Number of batches: ", x$nbatches, sep=""), "\n")
  cat(paste("Number(s) of observations (within each batch): ", paste(as.vector(table(x$batch)), collapse=", "), sep=""), "\n")
  cat(paste("Number of variables: ", ncol(x$xadj), sep=""), "\n")
  cat(paste("Number of factors: ", x$svobj$n.sv, sep=""), "\n")
  
}
