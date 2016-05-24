#' Standard covariance models for geostatistical data.
#' 
#' Creates a standard covariance model (\code{cmodStd}) object for geostatistical data.
#' 
#' The general form of the specified covariance function is \code{psill} * \eqn{\rho}(\code{d}; \code{r}) + (\code{evar} + \code{fvar})*(\code{d==0}), where \eqn{\rho} is the covariance function of the parametric models.
#' 
#' For the \code{exponential} model, \eqn{\rho}(\code{d}; \code{r}) is exp(-\code{d}/\code{r}).
#' 
#' For the \code{gaussian} model, \eqn{\rho}(\code{d}; \code{r}) is exp(-\code{d^2}/\code{r^2}).
#' 
#' For the \code{matern} model, \eqn{\rho}(\code{d}; \code{r}) is 2^(1-\code{par3})/\code{gamma}(\code{par3})*\code{sd}^\code{par3}*\code{besselK(sd, nu = par3)}, where \code{sd = d/r}.
#' 
#' For the \code{amatern} (alternative Matern) model, \eqn{\rho}(\code{d}; \code{r}) is \code{2^(1-par3)/gamma(par3)*sd^par3*besselK(sd, nu = par3)}, where \code{sd = 2 * sqrt(par3) * d/r}.
#' 
#' For the \code{spherical} model, \eqn{\rho}(\code{d}; \code{r}) is \code{1 - 1.5*sd + 0.5*(sd)^3} if \code{d < r}, and 0 otherwise, with \code{sd = d/r}.
#' 
#' For the \code{wendland1} model, \eqn{\rho}(\code{d}; \code{r}) is \code{(1 - sd)^4 * (4*sd + 1)} if \code{d < r}, and 0 otherwise, with \code{sd = d/r}.
#' 
#' For the \code{wendland2} model, \eqn{\rho}(\code{d}; \code{r}) is \code{(1 - sd)^6 * (35*sd^2 + 18*sd + 3))/3} if \code{d < r}, and 0 otherwise, with \code{sd = d/r}.
#' 
#' For the \code{wu1} model, \eqn{\rho}(\code{d}; \code{r}) is \code{(1 - sd)^3 * (1 + 3*sd + sd^2)} if \code{d < r}, and 0 otherwise, with \code{sd = d/r}.
#' 
#' For the \code{wu2} model, \eqn{\rho}(\code{d}; \code{r}) is \code{(1 - sd)^4*(4 + 16*sd + 12*sd^2 + 3*sd^3))/4} if \code{d < r}, and 0 otherwise, with \code{sd = d/r}.
#' 
#' For the \code{wu3} model, \eqn{\rho}(\code{d}; \code{r}) is \code{(1 - sd)^6 * (1 + 6*sd + 41/3*sd^2 + 12*sd^3 + 5*sd^4 + 5/6*sd^5)} if \code{d < r}, and 0 otherwise, with \code{sd = d/r}.
#' 
#' @param model A standard semivariance model type.
#' @param psill The partial sill of the model.  Must be a postive number.
#' @param r The range parameter r.  Must be a positive number.
#' @param evar The variance of the errors.  Must be non-negative number.  The default is 0.
#' @param fvar The finescale variance (microscale error).  Must be a non-negative number.  The default is 0.
#' @param par3 The value of the third parameter for 3 parameter models.  Must be a positive number.  The default is 0.5.
#' 
#' @return Returns a \code{cmodStd} object.
#' 
#' @author Joshua French
#' @export
#' @references Waller, L. A., & Gotway, C. A. (2004). Applied Spatial Statistics for Public Health Data. John Wiley & Sons.
#' @seealso \code{\link[spam]{covmat}}
#' @examples 
#' cmod.std(model = "exponential", psill = 1, r = 1) 

cmod.std = function(model, psill, r, evar = 0, 
                    fvar = 0, par3 = 0.5)
{
  check.args.cmod.std(model, psill, r, evar, fvar, par3)
  out = list(model = model, psill = psill, r = r, 
             evar = evar, fvar = fvar, par3 = par3)
  class(out) = "cmodStd"
  return(out)
}

check.args.cmod.std = function(model, psill, r, evar, 
                                fvar, par3)
{
  if(!is.numeric(psill) || 
     length(psill) != 1 || 
     min(psill) <= 0 || !is.finite(psill))
  {
    stop("psill must be a finite positive number.")
  }
  if(!is.numeric(r) || length(r) != 1 || 
     min(r) <= 0 || !is.finite(r))
  {
    stop("r must be a finite positive number.")
  }
  if(!is.numeric(evar) || length(evar) != 1 || 
     min(evar) < 0 || !is.finite(evar))
  {
    stop("evar must be a finite non-negative number.")
  }
  if(!is.numeric(fvar) || length(fvar) != 1 || 
     min(fvar) < 0 || !is.finite(fvar))
  {
    stop("fvar must be a non-negative number.")
  }
  if(!is.numeric(par3) || length(par3) != 1 || 
     min(par3) <= 0)
  {
    stop("par3 must be a finite positive number.")
  }
}

