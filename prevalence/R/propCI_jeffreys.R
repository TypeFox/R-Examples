propCI_jeffreys <-
function(x, n, l){

  lw <- qbeta(l[1], x + .5, n - x + .5)
  up <- qbeta(l[2], x + .5, n - x + .5)

  return(c(lw, up))
}