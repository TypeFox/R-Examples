#' Finds maximum likelihood estimates of model parameters for a geostatistical model
#' 
#' \code{mle} estimates the parameters of a geostatistical linear model.  The function is written to automatically adapt based on the class of \code{object}.  See Details.
#' 
#' In the case of a \code{geolmStd} \code{object}, the likelihood has been concentrated so that only the range parameter \code{r} and a scale parameter \code{lambda = nugget/psill} need to be estimated.
#' 
#' If \code{object} is a \code{geolmStd}, then \code{lower} is of length 2 if the covariance model of \code{cmod} is not \code{matern} or \code{amatern}.  Otherwise, it should be of length 3.  The first parameter is related to the range parameter \code{r}, the second to the scale parameter \code{lambda}, and the third to \code{par3}, if applicable.  If \code{lower = NULL}, then the lower bounds are 0.001, 0, and 0.1, respectively.  A similar pattern holds for \code{upper}, with the default being \code{3 * max(d)}, where \code{d} is the matrix of distances between coordinates, \code{5}, and \code{2.5}.
#' 
#' The \code{kkt} argument in the \code{control} list is set to be \code{FALSE}.
#' 
#' 
#' @param object A geostatistical linear model object producted by the \code{geolm} function.
#' @param reml  A logical value indicating whether standard maximum likelihood estimation should be performed (\code{reml = FALSE}).  If \code{reml = TRUE}, then restricted maximum likelihood is performed.  Defaul is \code{FALSE}.
#' @param est A character vector indicator whether the error variance (\code{est="e"}) or finescale variance (\code{est = "f"}) should be estimated.  The other component of the nugget variance is held constant, and in the case of a \code{geolmStd} object, is set to 0.
#' @param lower A vector of 2 or 3 specifying the lowerbound of parameter values.  See Details.
#' @param upper lower A vector of 2 or 3 specifying the lowerbound of parameter values.  See Details.
#' @param method The optimization method.  Default is \code{"nlminb"}, with \code{"L-BFGS-B"} being another acceptable choice.  See \code{\link[optimx]{optimx}} for details.
#' @param itnmax An integer indicating the maximum number of iterations to allow for the optimization prodedure.  
#' @param control A list of control parameters passed internally to \code{\link[optimx]{optimx}}.  See \code{\link[optimx]{optimx}} for details.
#' @param ... Currently unimplemented.
#' 
#' @author Joshua French
#' @export
#' @examples 
#' set.seed(10)
#' n = 100
#' @rdname mle
#' @export
mle = function(object, reml = FALSE, est = "e", ...)
{
  UseMethod("mle")
}

#' @rdname mle
#' @importFrom optimx optimx
#' @export
mle.geolmStd <- function(object, reml = FALSE, est = "e", 
                        lower = NULL, upper = NULL, method = "nlminb",
                        itnmax = NULL, control = list(), ...)
{
  if(!is.logical(reml)) stop("reml should be a logical value")
  if(!is.element(est, c("e", "f"))) stop("invalid choice for est")
  
  scmod = object$cmod
  scmod$psill = 1
  scmod$evar = 0
  matern = is.element(scmod$model, c("matern", "amatern"))
  if(est == "e")
  {
    weights = object$vediag/object$evar
    lambda = object$evar/object$cmod$psill
    scmod$fvar = 0
  }else
  {
    weights = rep(1, length(object$y))
    lambda = object$cmod$fvar/object$cmod$psill
    scmod$evar = 0
  }
  d = sp::spDists(as.matrix(object$coords), longlat = object$longlat)
  
  # if there are no lower bounds and the model is not a matern
  if(is.null(lower) & !is.element(scmod$model, c("matern", "amatern")))
  {
    lower = c(0.01, 0)
  }else
  {
    lower = c(0.01, 0, 0.1)
  }
  # same thing for upper bounds
  if(is.null(upper) & !is.element(scmod$model, c("matern", "amatern")))
  {
    upper = c(max(d) * 3, 5)
  }else
  {
    upper = c(max(d) * 3, 5, 2.5)
  }
  
  ini_par = c(object$cmod$r, lambda)
  if(is.element(scmod$model, c("matern", "amatern")))
  {
    ini_par = c(ini_par, object$cmod$par3)
  }
  control$kkt = FALSE
  
  out = optimx::optimx(ini_par, fn = ploglik.cmodStd, 
                 lower = lower, upper = upper,
                 y = object$y, x = object$x, d = d,
                 weights = weights,  
                 scmod = scmod, 
                 reml = reml,
                 return.ll = TRUE,
                 method = method,
                 itnmax = itnmax, control = control)

  ploglik.cmodStd(ini_par, y = object$y, x = object$x, d = d, 
                             weights = weights, scmod = scmod, reml = reml, 
                             return.ll = TRUE)
    
  
  if(length(ini_par) == 2)
  {
    new_par = c(out$p1, out$p2)
    psill = ploglik.cmodStd(new_par, y = object$y, x = object$x, weights = weights, scmod = scmod, d = d, reml = reml, return.ll = FALSE)
  }else
  {
    new_par = c(out$p1, out$p2, out$p3)
    psill = ploglik.cmodStd(new_par, y = object$y, x = object$x, weights = weights, scmod = scmod, d = d, reml = reml, return.ll = FALSE)
  }
  
  new_cmod = object$cmod_evar0
  
  # determine whether nugget is error or finescale, replace appropriately
  if(est == "e"){
    new_cmod$evar = new_par[2]*psill
  }else
  {
    new_cmod$fvar = new_par[2]*psill
  }
  new_cmod$psill = psill
  new_cmod$r = out$p1
  if(matern) new_cmod$par3 = out$p3
  
  newgeolm = update.geolmStd(object, new_cmod)
  newgeolm$optimx = out
  return(newgeolm)
}

# return -2 * log likelihood (if return.ll = TRUE), or the estimate
# of the psill, if return.ll = FALSE.
ploglik.cmodStd = function(parm, y, x, d, weights, scmod, reml, return.ll = TRUE)
{
  scmod$r = parm[1]
  if(length(parm) > 2) scmod$par3 = parm[3]
  v = eval.cmod(scmod, d) + parm[2] * diag(1/weights)
  n = nrow(v)
  cholv = chol(v)
  
  cholvix = backsolve(cholv, x, transpose = TRUE)
  vix = forwardsolve(cholv, cholvix, upper.tri = TRUE)
  xtvix = crossprod(cholvix)
  cholxtvix = chol(xtvix)
  map2coeff = forwardsolve(cholxtvix, 
                           backsolve(cholxtvix, t(vix), 
                                     transpose = TRUE), 
                           upper.tri = TRUE)
  coeff = map2coeff %*% y
  resid = y - x %*% coeff
  cholviresid = backsolve(cholv, resid, transpose = TRUE)
  sigmasq = crossprod(cholviresid)/n
  
  # mle likelihoods in long form
  # -1/2*(n * log(2*pi) + determinant(sigmasq[1,1] * v, log = TRUE)$mod + n)
  # -1/2*(n * log(2*pi) + determinant(sigmasq[1,1] * v, log = TRUE)$mod + 
  # crossprod(cholviresid)/sigmasq[1,1])
  
  # -1/2*(n * log(2*pi) + determinant(sigmasq[1,1] * v, log = TRUE)$mod + 
  #         determinant(xtvix/sigmasq[1, 1], log = TRUE)$mod + 
  #          crossprod(cholviresid)/sigmasq[1,1])
  
  # both are -2 the loglikelihood
  if(!reml)
  {
    ll = (n * log(2*pi) + n * log(sigmasq) + 
               sum(2 * log(diag(cholv))) + n)
  }else
  {
    sigmasq = sigmasq *n/(n - ncol(x))
    ll = (n * log(2*pi) + (n - ncol(x)) * log(sigmasq) + 
                        sum(2 * log(diag(cholv))) + 
                        determinant(xtvix, logarithm = TRUE)$mod
                        + n - ncol(x))
  }
  # reml likelihood in long form
  # sigmasq2 = crossprod(cholviresid)/(n - ncol(x))
  # -1/2*(n * log(2*pi) + determinant(sigmasq2[1,1] * v, log = TRUE)$mod + 
  #          determinant(xtvix/sigmasq2[1, 1], log = TRUE)$mod + 
  #        crossprod(cholviresid)/sigmasq2[1,1])
  if(return.ll){
    return(ll[1,1])
  }else
  {
    return(sigmasq[1, 1])
  }
}