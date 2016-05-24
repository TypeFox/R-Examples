rad2mat <-
function(pyDat, fill.N = TRUE) {
## shoves RAD sequence into an individuals x loci matrix
  if(!"radSummary" %in% names(pyDat)) pyDat$radSummary <- summary(pyDat) # calls summary if it wasn't done at read-time 
  loci <- dimnames(pyDat$radSummary$inds.mat)[[2]]
  inds <- dimnames(pyDat$radSummary$inds.mat)[[1]]
  out <- matrix(NA, length(inds), length(loci), dimnames = list(inds, loci))
  counter = 0
  for(i in loci) {
    if(counter %/% 1000 == counter / 1000) message(paste('... making rad.mat for', i))
	counter <- counter + 1
	out[, i] <- pyDat$radSummary$seqs.per.locus[[i]][inds]
    if(fill.N) out[inds[!inds %in% pyDat$radSummary$tips.per.locus[[i]]], i] <- paste(rep("N", pyDat$radSummary$locus.lengths[[i]]), collapse = "")
	}
  class(out) <- 'rad.mat'
  return(out)
}
