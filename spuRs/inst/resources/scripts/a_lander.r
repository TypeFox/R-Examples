a_lander <- function(amax) {
  # Find that thrust a and firing time Ti which land the rocket using
  # the least fuel. amax is the maximum thrust possible
  g <- function(a) {
    Ti <- Ti_lander(a)
    return(lander(Ti, a)[3])
  }
  opt <- optimize(g, c(0, amax), maximum=T)
  a <- opt$maximum     # optimal thrust
  mf <- opt$objective  # final mass
  return(c(Ti_lander(a), a, mf))
}