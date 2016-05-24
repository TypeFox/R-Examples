# ========================================================================
# .sims.as.list  -  split the simulations into a list
# ========================================================================

.sims.as.list <- function (x) {
  ## retain dimensions, and always return getnsims() simulations.
  if (is.null(x)) return(NULL)
  s <- sims(as.rv(x), n.sims=getnsims())
  s <- split(s, row(s)) ## faster than applying list to 1=rows.
  if (!is.null(d <- dim(x))) {
    dn <- dimnames(x)
    s <- lapply(s, function (x) { dim(x) <- d; dimnames(x) <- dn; x})
  }
  ## The default names will be "1", "2", ... etc.
  ## set names to NULL since this may interfere with mapply( ... ) in "[<-.rv"
  names(s) <- NULL
  return(s)
}
