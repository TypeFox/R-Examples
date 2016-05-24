
SampleStratified <- function(idxTrue, scale=1, verbose=TRUE) {
  #
  # given an index of true values, returns an index with stratified sampling
  #
  stopifnot(class(idxTrue) == "logical")
  nTrue  <- sum( idxTrue)
  nFalse <- sum(!idxTrue)
  if (verbose) {
    cat("Executing stratified sampling:\n")
    cat(sprintf("  Before: %d records, %d / %d true / false, %8.6f true rate\n", 
            length(idxTrue), nTrue, nFalse, nTrue/length(idxTrue) ))
  }
  sampleRate <- sqrt(nFalse / nTrue) / scale # sqrt ratio of false to true
  # if rate is < 1 then there are more true then false, return all rows
  if (sampleRate < 1) { 
    return(1:length(idxTrue)) 
  }
  # get indices of false rows
  idxFalseAll <- which(!idxTrue)
  numKeep <- round(nFalse / sampleRate) # number of false elements to keep
  idxFalseKeep <- sample( idxFalseAll, numKeep ) # a random sample of the false indices
  idxTrueKeep <- which(idxTrue)
  idxKeep <- append(idxFalseKeep, idxTrueKeep) # indices of rows to keep, unsorted
  idxKeep <- sort(idxKeep) # sorted ascending indexes
  if (verbose) {
    cat(sprintf("  After : %d records, %d / %d true / false, %8.6f true rate\n", 
            length(idxKeep), nTrue, length(idxFalseKeep), nTrue/length(idxKeep) ))
  }
  logicalKeep <- rep(FALSE, length(idxTrue))
  logicalKeep[idxKeep] <- TRUE
  return( logicalKeep ) # return logical vector of records to keep
}