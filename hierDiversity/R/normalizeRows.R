normalizeRows <-
function(dat) {
  if (is.matrix(dat)) {
    norm.abundance <- dat/rowSums(dat)
    norm.abundance[is.nan(norm.abundance)] <- 0
    norm.abundance
  }
  if (is.vector(dat)) {
    norm.abundance <- dat/sum(dat)
    norm.abundance[is.nan(norm.abundance)] <- 0
  }
  return(norm.abundance)
}
