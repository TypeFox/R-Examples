#' Predict method for geostatistical models
#' 
#' \code{predict} calculates the predicted values at specified locations.  The method can additionally provide the mean square prediction error (mspe) and perform conditional simulation.
#' 
#' The \code{newdata} data frame must include the relevant covariates for the prediction locations, where the covariates are specified on the right side of the \code{~} in \code{object$formula}.  \code{newdata} must also include the coordinates of the prediction locations, with these columns having the names provided in \code{object$coordnames}.
#' 
#' @param object An object produced by the \code{geolm} function.
#' @param newdata An optional data frame in which to look for the coordinates at which to predict. If omitted, the observed data locations are used.
#' @param nsim A non-negative integer indicating the number of realizations to sample at the specified coordinates using conditional simulation.
#' @param vop The cross-covariance matrix between the observed responses and the responses to predict.
#' @param vp The covariance matrix of the responses to predict.
#' @param sp A logical value indicating whether to object returned should be of class \code{\link[sp]{SpatialPointsDataFrame}} for easier plotting with the \code{sp} package.  Default is \code{TRUE}.
#' @param dmethod The method used to decompose the covariance matrix for conditional simulation.  Valid options are \code{"chol"}, \code{"eigen"}, and \code{"svd"}.  The default is \code{"chol"}.
#' @param ... Currently unimplemented.
#' 
#' @return If \code{sp = TRUE}, then a \code{SpatialPointDataFrame} from the \code{sp} package is returned, with components including the prediction coordinates, the predicted responses \code{pred}, the mean square prediction error (\code{mspe}), the root mean square prediction error (\code{rmspe}), and the conditional realizations, if application (\code{sim.1}, \code{sim.2}, ...).  If \code{sp = FALSE}, then a list of class \code{gearKrige} is returned, with components \code{pred}, \code{mspe}, \code{rmspe}, and \code{sim}, if relevant.   
#' @author Joshua French
#' @importFrom stats predict delete.response rnorm model.matrix terms
#' @examples 
#' 
#' # generate response
#' y = rnorm(10)
#' # generate coordinates
#' x1 = runif(10); x2 = runif(10)
#' 
#' # data frame for observed data
#' data = data.frame(y, x1, x2)
#' coords = cbind(x1, x2)
#' d = as.matrix(dist(coords))
#' psill = 2 # partial sill
#' r = 4 # range parameter
#' evar = .1 # error variance
#' fvar = .1 # add finescale variance
#' # one can't generally distinguish between evar and fvar, but
#' # this is done for illustration purposes
#' 
#' # manually specify a an expoential covariance model 
#' v = psill * exp(-d/r) + (evar + fvar) * diag(10)
#' 
#' cmod_man = cmod.man(v = v, evar = evar)
#' 
#' #' # geolm for universal kriging
#' gearmod_uk = geolm(y ~ x1 + x2, data = data,
#'                  coordnames = c("x1", "x2"),
#'                  cmod = cmod_man)
#'
#' # newdata must have columns with prediction coordinates
#' # add 5 unsampled sites to sampled sites
#' newdata = data.frame(x1 = c(x1, runif(5)), x2 = c(x2, runif(5)))
#' newcoords = newdata[,c("x1", "x2")]
#' # create vop and vp using distances
#' dop = sp::spDists(as.matrix(coords), as.matrix(newcoords))
#' dp = as.matrix(dist(newcoords))
#' 
#' vop = psill * exp(-dop/r) + fvar * (dop == 0)
#' vp = psill * exp(-dp/r) + fvar * diag(nrow(newcoords))
#' 
#' # prediction for universal kriging, with conditional simulation, 
#' # using manual covariance matrices
#' pred_uk_man = predict(gearmod_uk, newdata, nsim = 2, vop = vop, vp = vp,
#'                       dmethod = "svd")
#' 
#' # do the same thing, but using the std covariance function
#' 
#' # prediction for universal kriging, with conditional simulation
#' cmod_std = cmod.std("exponential", psill = psill, r = r, evar = evar, fvar = fvar)
#' gearmod_uk2 = geolm(y ~ x1 + x2, data = data, coordnames = c("x1", "x2"), 
#'                     cmod = cmod_std)
#' pred_uk_std = predict(gearmod_uk2, newdata, nsim = 2, dmethod = "svd")
#' 
#' # compare results
#' range(pred_uk_man$pred - pred_uk_std$pred)
#' range(pred_uk_man$mspe - pred_uk_std$mspe)
#' @rdname predict.geolmMan
#' @export
predict.geolmMan = function(object, newdata, nsim = 0, vop,
                            vp, sp = TRUE, dmethod = "chol", ...) 
{
  #check validity of arguments
  if(nsim < 0 || !is.finite(nsim))
  { stop("nsim should be a non-negative integer") }
  if(nrow(vop) != length(object$y)) stop("nrow(vop) != length(object$y)")
  if(nrow(vp) != ncol(vp)) stop("vp must be a square matrix")
  if(ncol(vop) != nrow(vp)) stop("nrow(vop) != nrow(vp)")
  if(!is.element(dmethod, c("chol", "eigen", "svd")))
  {
    stop("Invalid dmethod type")
  }

  newcoords = as.matrix(newdata[,object$coordnames])

  # unneed if simple kriging
  if(is.null(object$mu))
  {
    newf = stats::delete.response(stats::terms(object$formula))
    newx = stats::model.matrix(newf, data = newdata)
    dimnames(newx) = NULL
  }
  
  # generate simulated data at observed and prediction locations
  if(nsim > 0)
  {
    n = nrow(object$coords); m = nrow(newcoords)
    
    # slightly different algorithm for simple kriging
    mu = 0
    if(!is.null(object$mu)) mu = object$mu
    
    newsim = decomp.cov(rbind(
      cbind(object$v - diag(object$vediag),vop),
      cbind(t(vop), vp)), method=dmethod) %*% 
      matrix(rnorm((n+m)*nsim), ncol = nsim) + mu
    
    # update various object to include simulated data
    object$y = cbind(object$y, newsim[1:n,] + 
                matrix(stats::rnorm(n*nsim, 
                       sd = sqrt(object$vediag)),
                       nrow = n, ncol = nsim))
    
    if(!is.null(object$mu))
    {
      object$cholviresid = 
        backsolve(object$cholv, object$y - mu, 
                  transpose = TRUE)
    }else
    {
      object$coeff = object$map2coeff %*% object$y
      object$cholviresid = backsolve(object$cholv, 
                                     object$y - object$x %*% object$coeff,
                                     transpose = TRUE)
    }
  }
  
  cholvivop = forwardsolve(object$cholv, vop, 
                           transpose = TRUE, 
                           upper.tri = TRUE)
  
  # simple kriging
  if(!is.null(object$mu))
  {
    pred = object$mu + crossprod(cholvivop, object$cholviresid)
    mspe = diag(vp) - colSums(cholvivop^2)
    mspe[mspe < 0] = 0 # numerical imprecision correction
  }else
  {
    pred = c(newx %*% object$coeff) + crossprod(cholvivop, object$cholviresid)

    mspe = diag(vp) - 
            colSums(cholvivop^2) + 
            colSums(forwardsolve(object$cholxtvix, 
                    t(newx - crossprod(cholvivop, 
                                       object$cholvix)), 
                    transpose = TRUE, upper.tri = TRUE)^2)
    mspe[mspe < 0] = 0 # numerical imprecision correction
  }
  
  # return results
  if(sp) # if return a SpatialPointsDataFrame
  {
    if(nsim > 0)
    {
      return(sp::SpatialPointsDataFrame(
        coords = newcoords, 
        data = data.frame(pred = pred[,1], 
                          mspe = mspe, 
                          rmspe = sqrt(mspe),
                          sim = newsim[-(1:n),] + 
                            (pred[,1] - pred[,-1]))))
    }else
    {
      return(sp::SpatialPointsDataFrame(
        coords = newcoords, 
        data = data.frame(pred = pred[,1], 
                          mspe = mspe, 
                          rmspe = sqrt(mspe))))
    }
  }else # if not SpatialPointsDataFrame
  {
    if(nsim > 0)
    {
      return(structure(
      list(pred = pred[,1], 
           mspe = mspe, 
           rmspe = sqrt(mspe),
           sim = newsim[-(1:n),] + (pred[,1] - pred[,-1])), 
      class = "gearKrige"))
    }else
    {
      return(structure(
        list(pred = pred[,1], 
             mspe = mspe, 
             rmspe = sqrt(mspe)), 
        class = "gearKrige"))
    }
  }
}

