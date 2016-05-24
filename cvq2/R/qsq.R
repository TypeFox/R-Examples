qsq <-
function( modelData, predictData, formula = NULL, nu = 0, round = 4, extOut = FALSE, extOutFile = NULL ){
#  call <- match.call()
  N <- nrow(predictData)

  #forward to the actual method
  x <- mainfunc.q2( modelData, predictData, formula, N, 1, nu, round, extOut, extOutFile, match.call() )
  return(x)
}

