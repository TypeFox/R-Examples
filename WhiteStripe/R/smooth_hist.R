#' @title Histogram smoothing for whitestripe
#'
#' @description Uses a generalized additive model (GAM) to smooth a 
#' histogram for whitestripe
#' @param x values of midpoints from \code{\link{hist}}
#' @param y values of counts from \code{\link{hist}}
#' @param deg degree of polynomials used
#' @param k Number of knots
#' @param method Method for smoothing for GAM
#' @param ... Arguments passed to \code{\link{gam}}
#' @export
#' @seealso \link[mgcv]{gam}
#' @importFrom mgcv gam
#' @return List of objects: x and y coordinates of histogram, coefficients from GAM, 
#' fitted values from GAM, the GAM model, the knots fittted, and degrees of polynomials
#' @examples 
#' data(t2.voi.hist)
#' y = t2.voi.hist$counts
#' x = t2.voi.hist$mids
#' x = x[!is.na(y)];
#' y = y[!is.na(y)]
#' # 30 used for speed of example
#' s.hist = smooth_hist(x, y, k=30)
#' plot(t2.voi.hist, border="red")
#' lines(s.hist)
smooth_hist = function(x, y, 
                       deg = 4, 
                       k = floor(min(250,length(x)/2)), 
                       method = "REML", ...){
  
  stopifnot(deg > 0)
#   print(paste0("Number of Knots: ", k))
  keep = which(y != -Inf)
  x = x[keep]
  y = y[keep]  
  
  qtiles <- seq(0, 1, length = k+2)[-c(1, k+2)]
  knots <- quantile(x, qtiles)
  phi = cbind(sapply(0:deg, function(k) (x^k) ), 
    sapply(knots, function(k) ((x - k > 0) * (x - k)^deg))
    )
  
  X = phi
  D = list(length=1)
  D[[1]] = diag(c(rep(0, deg+1), rep(1, k)))
  
#   print("Dimensions of D:", dim(D[[1]]))
  fit = gam(y~X-1, paraPen=list(X=D), method = method, ...)
  
  coefs = fit$coef
  fitted.vals = cbind(phi)%*%coefs
  
  ret = list(x, y, coefs, fitted.vals, fit, knots, deg)
  names(ret) = c("x", "y", "coefs", "fitted.vals", "fit", "knots", "deg")
  return(ret)
}



#' @title Gets $n^{th}$ derivative of smoothed histogram
#'
#' @description This function outputs the nth derivative of a histogram smooth.
#' @param x values from smooth_hist
#' @param coefs Coefficients from GAM from smooth_hist
#' @param knots Number of knots fit for GAM
#' @param deg Degree of polynomials
#' @param deriv.deg <what param does>
#' @export
#' @return Derivative of smoothed histogram
#' @examples 
#' data(smoothed_histogram)
#' dy<-get.deriv.smooth.hist(xvals, 
#' coefs=s.hist$coefs,
#' knots=s.hist$knots,
#' deg=s.hist$deg,
#' deriv.deg=1)
get.deriv.smooth.hist <- function(x,
                                  coefs,
                                  knots,
                                  deg=4,
                                  deriv.deg=1) {
  stopifnot(deg > 0)
  stopifnot(deriv.deg > 0)
  deriv.coefs<-coefs[2:length(coefs)]*c(1:deg,rep(deg,length(knots)))
  deriv.phi = cbind(
    sapply(0:(deg-1), function(k) (x^k) ), 
    sapply(knots, function(k) ((x - k > 0) * (x - k)^(deg-1)))
  )
  if (deriv.deg>1) {
    return(get.deriv.smooth.hist(x,coefs=deriv.coefs,deg=deg-1,knots=knots,deriv.deg=deriv.deg-1))
  } else {
    return(cbind(deriv.phi)%*%deriv.coefs)
  }
}