`motif.find` <-
function(motif, sequence) {
  ## return indices of motif within sequence
  position <- regexpr( paste(motif, collapse=""), paste(sequence,collapse=""))
  inds <- c(position):c(position+attr(position, "match.length")-1)
  return(inds)
}

