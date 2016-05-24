whiteXwhiteKernGradient <-
function (whiteKern1, whiteKern2, x1, x2, covGrad) {
  if ( nargs()<5 ) {
    covGrad <- x2
    x2 <- x1
  }
  
  if ( dim(as.matrix(x1))[2]>1 | dim(as.matrix(x2))[2]>1 )
    stop("Input can only have one column.")

  g <- list(g1=0, g2=0)
  return (g)
}
