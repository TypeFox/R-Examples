
##==============================================================================
## The Fiadeiro and Veronis scheme: advective finite difference weights (AFDW)
##==============================================================================

fiadeiro <- function(v, D, dx.aux=NULL, grid=list(dx.aux=dx.aux)) {

  Pe    <- abs(v)*grid$dx.aux/D
  sigma <- (1+(1/tanh(Pe)-1/Pe))/2
  return(sigma)
}

