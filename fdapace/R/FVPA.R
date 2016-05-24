#' Functional Variance Process Analysis for sparse or dense functional data
#' 
#' @param y A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param q A scalar defining the percentile of the pooled sample residual sample used for adjustment before taking log (default: 0.1).
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details in ?FPCA'.
#'
#'
#' @return A list containing the following fields:
#' \item{sigma2}{Variance estimator of the functional variance process.} 
#' \item{fpcaObjY}{FPCA object for the original data.} 
#' \item{fpcaObjR}{FPCA object for the functional variance process associated with the original data.} 
#' 
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10) 
#' #fvpaObj <- FVPA(sampWiener$yList, sampWiener$tList)
#' @references
#' \cite{Hans-Georg Mueller, Ulrich Stadtmuller and Fang Yao, "Functional variance processes." Journal of the American Statistical Association 101 (2006): 1007-1018}
#' 

FVPA = function(y, t, q= 0.1, optns = list()){ 
  
  if( (q <0) || (1 < q) ){
    warning("The value of 'q' is outside [0,1]; reseting to 0.1.")
  }
  
  fpcaObjY <- FPCA(y, t, optns)
  
  if( fpcaObjY$optns$dataType == 'Dense' ){
    stop("The data are considered to be 'dense' and therefore noiseless; FVPA is irrelevant.")
  }
  
  yFitted <- fitted(fpcaObjY);
  r = GetVarianceProcess(y, t, yFitted, workGrid = fpcaObjY$workGrid, delta = 0, logarithm = FALSE )
  delta = quantile(unlist(r), q)
  r = GetVarianceProcess(y, t, yFitted, workGrid = fpcaObjY$workGrid, delta = delta, logarithm = TRUE )
  
  fpcaObjR = FPCA(r, t, optns);
  return( list( sigma2 = fpcaObjR$sigma2, fpcaObjY = fpcaObjY, fpcaObjR = fpcaObjR))
}

GetVarianceProcess = function(y, t, yFitted, workGrid, delta =0.0, logarithm = FALSE){
  r <- list()
  for (i in 1:nrow(yFitted)){
    tempV = delta + (  y[[i]] - approx(x = workGrid, y = yFitted[i,], xout = t[[i]])$y )^2
      if(logarithm){
        r[[i]] = log(tempV);
      } else {
        r[[i]] = tempV;
      }
  }
  return(r)
}
