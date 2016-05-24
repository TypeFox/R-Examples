## Calculates the theoretical SSCM
# input: shape matrix (symmetric and pxp)
# output:
# theoretical SSCM


Shape2SSCM <- function(V) {
  p <- ncol(V)
  EVD <- eigen(V,symmetric=TRUE)
  delta <- evShape2evSSCM(EVD$values)
  sscm <- EVD$vectors %*% diag(delta) %*% t(EVD$vectors)
  return(sscm)
}
