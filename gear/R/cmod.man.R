#' Standard covariance models for geostatistical data.
#' 
#' \code{cmod.man} manually creates a covariance matrix object (\code{cmodMan}) for geostatistical data.
#' 
#' Note that \code{v} includes the error variance, i.e., \code{v = vz + ve}, where \code{vz} is the covariance matrix of the filtered (non-noisy) process, and the variance matrix of the errors is \code{ve = diag(evar/weights)}, where the weights come from the \code{geolm} object the covariance object is associated with.
#' 
#' @param v The covariance matrix of the observed data, including any errors.  The matrix should be square, symmetric, and positive definite, though that latter two conditions are not checked.
#' @param evar The variance of the errors.  Must be non-negative number.  The default is 0.
#' 
#' @return Returns a \code{cmodMan} object.
#' 
#' @author Joshua French
#' @export
#' @examples 
#' coords = matrix(runif(20), ncol = 2)
#' d = as.matrix(dist(coords))
#' cmod.man(v = exp(-d), evar = 1) 

cmod.man = function(v, evar = 0)
{
  check.args.cmod.man(v, evar)
  return(structure(list(v = v, evar = evar), class = "cmodMan"))
}

check.args.cmod.man = function(v, evar)
{
  if (!is.numeric(v) | nrow(v) != ncol(v) ) stop("v should be a square matrix")
  if (!is.numeric(evar) || length(evar) != 1 || 
     min(evar) < 0 || !is.finite(evar))
  {
    stop ("evar must be a finite non-negative number.")
  }
}

