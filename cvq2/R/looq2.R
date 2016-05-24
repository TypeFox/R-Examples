looq2 <-
function( modelData, formula = NULL, nu = 1, round = 4, extOut = FALSE, extOutFile = NULL  ){
#  call <- match.call()
  N <- nrow(modelData)

  #forward to the actual method
  x <- mainfunc.q2( modelData, NULL, formula, N, 1, nu, round, extOut, extOutFile, match.call() )
  return(x)
}

