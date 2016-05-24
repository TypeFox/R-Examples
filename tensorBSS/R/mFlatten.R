mFlatten <- function(x, m){
  rplusone <- length(dim(x))
  aperm(apply(x, c(m, rplusone), c), c(2, 1, 3))
}