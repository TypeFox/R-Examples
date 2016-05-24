`vmat` <-
function(wgths) {
  v <- as.matrix(wgths)
  r <-rowSums(v)              #row margins of weight matrix
  return(diag(r)-v)
}

