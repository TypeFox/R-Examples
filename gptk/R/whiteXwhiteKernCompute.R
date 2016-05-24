whiteXwhiteKernCompute <-
function (whiteKern1, whiteKern2, x1, x2) {
  if ( nargs()<4 )
    x2=x1

  k <- matrix(0, nrow=dim(as.array(x1))[1], ncol=dim(as.array(x2))[1])
  return (k)  
}
