locus.dist <-
function(pyIn, proportional = TRUE, upper = TRUE, diagonal = TRUE) {
  if(class(pyIn) == 'pyRAD.loci') pyIn <- pyIn$radSummary$inds.mat
  numInds <- dim(pyIn)[1]
  numLoci <- ifelse(proportional, dim(pyIn)[2], 1)
  out <- matrix(NA, numInds, numInds, dimnames = list(row.names(pyIn), row.names(pyIn)))
  for(i in 1:numInds) {
    for(j in 1:i) {
	  out[i, j] <- sum(colSums(pyIn[c(i, j), ]) == 2) / numLoci
	  }}
  if(upper) out <- as.matrix(as.dist(out))
  if(diagonal) diag(out) <- apply(pyIn, 1, sum) / numLoci
  class(out) <- 'locus.dist'
  out
  }
