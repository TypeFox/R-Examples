#' Builds response Y and predictor X.
#' 
#' This function builds the response variables Y and predictor variables X from
#' the input data.
#' 
#' 
#' @param targetData Target input data.
#' @param predData Predictor input data.
#' @param GLOBvar Global variables of the MCMC simulation.
#' @return A list with elements: \item{X}{Predictor variables.}
#' \item{Y}{Response variables.}
#' @author Sophie Lebre
#' @export buildXY
buildXY <-
function(targetData, predData, GLOBvar){
  ### Build response Y and predictor X

  ### Assignement of global variables used here ###
  n = GLOBvar$n
  m = GLOBvar$m
  q = GLOBvar$q
  dyn = GLOBvar$dyn
  target = GLOBvar$target
  bestPosMat = GLOBvar$bestPosMat
  ### end assignement ###

  # Read target data
  Y = as.array(readDataTS(data=targetData, posI=target, t0=dyn, tf=n, m=m, n=n))

  # Read predictor data
  posTF = as.matrix(bestPosMat)[target,1:q]
  dataTF = t(readDataTS(data=predData, posI=posTF, t0=0, tf=n-dyn, m=m, n=n))

  ## Add a constant vector to the predictor data
  X = cbind(dataTF,array(1,length(dataTF[,1])))
  
  # Return formatted data
  return(list(X=X,Y=Y))
}

