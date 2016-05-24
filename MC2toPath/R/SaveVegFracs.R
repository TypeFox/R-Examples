SaveVegFracs = function(vegChanges, outFile) {
  cat(c("annual fraction of cells in each vegetation type from file:", vegChanges[[1]], "\n"), file=outFile, append=FALSE)
  startYear = vegChanges[[2]][1]
  vts = vegChanges[[3]]
  cat(c("year (down) X VTYPE (across)"), file=outFile, append=TRUE)
  localMatrix = vegChanges[[4]]  
  nSeq = dim(localMatrix)[2] # length(sequence)
  stopifnot(nSeq>=1)
  
  matrixRows = dim(localMatrix)[1]
  stopifnot(matrixRows>1)
  
  for (ndx in 1:matrixRows) { cat(c(", ", vts[ndx]), file=outFile, append=TRUE) }
  cat(c("\n"), file=outFile, append=TRUE)
  
  indexedSeq = matrix(0, nrow=nSeq, ncol=matrixRows+1)
  
  for (ndx in 1:nSeq) {
    cat(c(startYear + ndx - 1), file=outFile, append=TRUE)
    indexedSeq[ndx, 1] = startYear+ndx-1
    for (row in 1:matrixRows) {
      if (is.na(localMatrix[row, ndx])) localMatrix[row, ndx] = 0
      cat(c(", ", localMatrix[row, ndx]), file=outFile, append=TRUE)
      indexedSeq[ndx, row+1] = localMatrix[row, ndx]
    } # end of for (row in 1:matrixRows)
    cat(c("\n"), file=outFile, append=TRUE)    
  } # end of for (ndx in 1:nSeq)
      
  return(TRUE)
}