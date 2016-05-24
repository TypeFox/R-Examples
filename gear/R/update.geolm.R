#' Update linear model for geostatistical data.
#' 
#' \code{update} updates a geostatistical linear model based on the given covariance model.
#' 
#' @param object An object produced by the \code{geolm} function.
#' @param cmod A covariance model object obtained from one of the \code{cmod.*} functions or the \code{mle} function.
#' @param ... Not implemented.
#' @return Returns an object of the same class as \code{object}.
#' 
#' @author Joshua French
#' @importFrom stats update
#' @importFrom sp spDists
#' @examples 
#' # generate response
#' y = rnorm(10)
#' # generate coordinates
#' x1 = runif(10); x2 = runif(10)
#' 
#' # data frame for observed data
#' data = data.frame(y, x1, x2)
#' coords = cbind(x1, x2)
#' psill = 2 # partial sill
#' r = 4 # range parameter
#' evar = .1 # error variance
#' fvar = .1 # add finescale variance
#' # one can't generally distinguish between evar and fvar, but
#' # this is done for illustration purposes
#' 
#' cmod_std = cmod.std("exponential", psill = psill, r = r, 
#'                     evar = evar, fvar = fvar)
#' 
#' cmod_std2 = cmod.std("exponential", psill = psill + 1, r = r + .5, 
#'                      evar = evar + .01, fvar = fvar)
#' 
#' # check geolm update for universal kriging
#' gear1 = geolm(y ~ x1 + x2, data = data,
#'               coordnames = c("x1", "x2"),
#'               cmod = cmod_std)
#' gear2 = geolm(y ~ x1 + x2, data = data,
#'               coordnames = c("x1", "x2"),
#'               cmod = cmod_std2)
#' gear2b = update(gear1, cmod_std2)
#' identical(gear2, gear2b)
#' @rdname update
#' @export
update.geolmStd = function(object, cmod, ...)
{
  if(class(cmod) != "cmodStd") stop("The class of cmod doesn't match the class of object")
  
  vediag = cmod$evar/object$weights
  # create covariance matrix for observed data

  cmod_evar0 = cmod
  cmod_evar0$evar = 0
    
  v = eval.cmod(cmod_evar0, 
              d = sp::spDists(object$coords, longlat = object$longlat)) + 
              diag(vediag)
  cholv = chol(v)
      
  if(is.null(object$mu))
  {
    ###compute matrix products for future use
    cholvix = backsolve(cholv, object$x, transpose = TRUE)
    vix = forwardsolve(cholv, cholvix, upper.tri = TRUE)
    xtvix = crossprod(cholvix)
    cholxtvix = chol(xtvix)
    vcov = chol2inv(cholxtvix)
    map2coeff = forwardsolve(cholxtvix, 
                             backsolve(cholxtvix, t(vix), 
                                       transpose = TRUE), 
                             upper.tri = TRUE)
    
    coeff = map2coeff %*% object$y
    resid = object$y - object$x %*% coeff
  }else
  {
    vix = NULL
    xtvix = NULL
    resid = object$resid
    coeff = NULL
    cholvix = NULL
    cholxtvix = cholxtvix = NULL 
    vcov = NULL
    map2coeff = NULL
  }
      
  cholviresid = backsolve(cholv, resid, transpose = TRUE)
  n = nrow(v)
  loglik = -n/2*log(2*pi) - 
    sum(2 * log(diag(cholv)))/2 - 
    crossprod(cholviresid)[1,1]/2

  return(structure(list(y = object$y, 
                        x = object$x, 
                        loglik = loglik,
                        resid = resid, 
                        coeff = coeff, 
                        mu = object$mu, 
                        v = v, 
                        cholv = cholv,
                        cholvix = cholvix, 
                        cholviresid = cholviresid,
                        cholxtvix = cholxtvix, 
                        vcov = vcov,
                        map2coeff = map2coeff, 
                        formula = object$formula, 
                        coordnames = object$coordnames,
                        cmod_evar0 = cmod_evar0,
                        weights = object$weights, 
                        vediag = vediag,
                        evar = cmod$evar,
                        coords = object$coords,
                        longlat = object$longlat), class = "geolmStd"))
}

#' @rdname update
#' @export
update.geolmMan = function(object, cmod, ...)
{
  if(class(cmod) != "cmodMan") stop("The class of cmod doesn't match the class of object")
  
  vediag = cmod$evar/object$weights
  # create covariance matrix for observed data
  
  cmod_evar0 = NULL
  v = cmod$v
  cholv = chol(v)
  
  if(is.null(object$mu))
  {
    ###compute matrix products for future use
    cholvix = backsolve(cholv, object$x, transpose = TRUE)
    vix = forwardsolve(cholv, cholvix, upper.tri = TRUE)
    xtvix = crossprod(cholvix)
    cholxtvix = chol(xtvix)
    vcov = chol2inv(cholxtvix)
    map2coeff = forwardsolve(cholxtvix, 
                             backsolve(cholxtvix, t(vix), 
                                       transpose = TRUE), 
                             upper.tri = TRUE)
    
    coeff = map2coeff %*% object$y
    resid = object$y - object$x %*% coeff
  }else
  {
    vix = NULL
    xtvix = NULL
    resid = object$resid
    coeff = NULL
    cholvix = NULL
    cholxtvix = cholxtvix = NULL 
    vcov = NULL
    map2coeff = NULL
  }
  
  cholviresid = backsolve(cholv, resid, transpose = TRUE)
  n = nrow(v)
  loglik = -n/2*log(2*pi) - 
    sum(2 * log(diag(cholv)))/2 - 
    crossprod(cholviresid)[1,1]/2
  
  return(structure(list(y = object$y, 
                        x = object$x, 
                        loglik = loglik,
                        resid = resid, 
                        coeff = coeff, 
                        mu = object$mu, 
                        v = v, 
                        cholv = cholv,
                        cholvix = cholvix, 
                        cholviresid = cholviresid,
                        cholxtvix = cholxtvix, 
                        vcov = vcov,
                        map2coeff = map2coeff, 
                        formula = object$formula, 
                        coordnames = object$coordnames,
                        cmod_evar0 = cmod_evar0,
                        weights = object$weights, 
                        vediag = vediag,
                        evar = cmod$evar,
                        coords = object$coords,
                        longlat = object$longlat), class = "geolmMan"))
}

