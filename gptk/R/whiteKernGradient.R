whiteKernGradient <-
function (kern, x, x2, covGrad) {
  if ( nargs()==3 ) {
    covGrad <- x2
    g <- sum(diag(as.matrix(covGrad)))
  } else {
    g <- 0
  }  

  return (g)
}
