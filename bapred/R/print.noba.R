print.noba <-
function(x, ...) {

  cat("Unadjusted training data.", "\n")
  cat(paste("Number of batches: ", x$nbatches, sep=""), "\n")
  cat(paste("Number(s) of observations (within each batch): ", paste(as.vector(table(x$batch)), collapse=", "), sep=""), "\n")
  cat(paste("Number of variables: ", ncol(x$xadj), sep=""), "\n")
  
}
