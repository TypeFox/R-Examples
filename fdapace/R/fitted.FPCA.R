#' Fitted functional sample from FPCA object
#' 
#' Combine the zero-meaned fitted values and the interpolated mean to get the fitted values for the trajectories or the derivatives of these trajectories.
#' 
#' @param object A object of class FPCA returned by the function FPCA().   
#' @param k The integer number of the first k components used for the representation. (default: length(fpcaObj$lambda ))
#' @param derOptns A list of options to control the derivation parameters specified by \code{list(name=value)}. See `Details'. (default = NULL)
#'
#' @details Available derivation control options are 
#' \describe{
#' \item{p}{The order of the derivatives returned (default: 0, max: 2)}
#' \item{method}{The method used to produce the sample of derivatives ('EIG' (default) or 'QUO'). See Liu and Mueller (2009) for more details}
#' \item{bw}{Bandwidth for smoothing the derivatives (default: p * 0.10 * S)}
#' \item{kernelType}{Smoothing kernel choice; same available types are FPCA(). default('epan')}
#' }
#' @param ... Additional arguments
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' fittedY <- fitted(res)
#' @references
#' \cite{Liu, Bitao, and Hans-Georg Mueller. "Estimating derivatives for samples of sparsely observed functions, with application to online auction dynamics." Journal of the American Statistical Association 104, no. 486 (2009): 704-717. (Sparse data FPCA)}
#' @export


fitted.FPCA <-  function (object, k = NULL, derOptns = list(), ...) {
  
  derOptns <- SetDerOptions(fpcaObject = object, derOptns)
  p <- derOptns[['p']]
  method <- derOptns[['method']]
  bw <-  derOptns[['bw']] # 
  kernelType <- derOptns[['kernelType']]

  fpcaObj <- object
  if (class(fpcaObj) != 'FPCA'){
    stop("fitted.FPCA() requires an FPCA class object as basic input")
  }

  if( is.null(k) ){
    k = length( fpcaObj$lambda )
  } else {
    if( ( round(k)>=1) && ( round(k) <= length( fpcaObj$lambda ) ) ){
      k = round(k);
    } else {
      stop("'fitted.FPCA()' is requested to use more components than it currently has available. (or 'k' is smaller than 1)")
    }
  }
 
  if( ! (p %in% c(0,1,2))){
    stop("'fitted.FPCA()' is requested to use a derivative order other than p = {0,1,2}!")
  } 

  if( p < 1 ){  
    ZMFV = fpcaObj$xiEst[,1:k, drop = FALSE] %*% t(fpcaObj$phi[,1:k, drop = FALSE]);   
    IM = fpcaObj$mu 
    return( t(apply( ZMFV, 1, function(x) x + IM))) 
  } else { #Derivative is not zero
 
    if( k > SelectK( fpcaObj, FVEthreshold=0.95, criterion='FVE')$k ){
    warning("Potentially you use too many components to estimate derivatives. \n  Consider using SelectK() to find a more informed estimate for 'k'.");
    }

    if( is.null(method) ){
      method = 'EIG'
    }

    mu = fpcaObj$mu
    phi = fpcaObj$phi
    obsGrid = fpcaObj$obsGrid
    workGrid = fpcaObj$workGrid

    if ( method == 'EIG'){
      phi = apply(phi, 2, function(phiI) Lwls1D(bw = bw, kernelType, win = rep(1, length(workGrid)), 
                                                  xin = workGrid, yin = phiI, xout = workGrid, npoly = p, nder = p))
      mu = Lwls1D(bw = bw, kernelType, win = rep(1, length(workGrid)), xin = workGrid, yin = mu, xout = workGrid, npoly = p, nder = p)
      ZMFV = fpcaObj$xiEst[,1:k, drop = FALSE] %*% t(phi[,1:k, drop = FALSE]);
      IM = mu 
      return( t(apply( ZMFV, 1, function(x) x + IM) ))
    }

    if( method == 'QUO'){
      impSample <- fitted(fpcaObj, k = k); # Look ma! I do recursion!
      return( t(apply(impSample, 1, function(curve) Lwls1D(bw = bw, kernelType, win = rep(1, length(workGrid)), 
                                                         xin = workGrid, yin = curve, xout = workGrid, npoly = p, nder = p))))
    }
  }

  stop('You asked for a derivation scheme that is not implemented.')
  return(NULL)
}

getEnlargedGrid <- function(x){
  N <- length(x)
  return (  c( x[1] - 0.1 * diff(x[1:2]), x, x[N] + 0.1 * diff(x[(N-1):N])) )
}

getDerivative <- function(y, t, ord=1){  # Consider using the smoother to get the derivatives
  if( length(y) != length(t) ){
    stop("getDerivative y/t lengths are unequal.")
  }
  newt = getEnlargedGrid(t) # This is a trick to get first derivatives everywhere
  newy = Hmisc::approxExtrap(x=t, y=y, xout= newt)$y

  if (ord == 1) {
    der <- numDeriv::grad( stats::splinefun(newt, newy) , x = t )
  } else if (ord == 2) {
    der <- sapply(t, function(t0) 
                  numDeriv::hessian( stats::splinefun(newt, newy) , x = t0 )
                  )
  }

  return(der)
}

getSmoothCurve <- function(t, ft, GCV = FALSE, kernelType = 'epan', mult = 1){
  myBw = ifelse( GCV, GCVLwls1D1( yy= ft, tt =t, npoly=1, nder=0, dataType='Sparse', kernel=kernelType)[['bOpt']] ,
                      CVLwls1D(   y= ft, t = t, npoly=1, nder=0, dataType='Sparse', kernel=kernelType, kFolds = 10))
  myBw <- myBw * mult
  smoothCurve = Lwls1D(bw = myBw, kernel_type= kernelType, win = rep(1, length(t)), yin = ft, xout = t, xin= t)
  return(smoothCurve)
}




