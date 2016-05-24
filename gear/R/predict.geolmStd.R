#' Predict method for geostatistical models
#' 
#' \code{predict} calculates the predicted values at specified locations.  The method can additionally provide the mean square prediction error (mspe) and perform conditional simulation.
#' 
#' The \code{newdata} data frame must include the relevant covariates for the prediction locations, where the covariates are specified on the right side of the \code{~} in \code{object$formula}.  \code{newdata} must also include the coordinates of the prediction locations, with these columns having the names provided in \code{object$coordnames}.
#' 
#' @param object An object produced by the \code{geolm} function.
#' @param newdata An optional data frame in which to look for the coordinates at which to predict. If omitted, the observed data locations are used.
#' @param nsim A non-negative integer indicating the number of realizations to sample at the specified coordinates using conditional simulation.
#' @param sp A logical value indicating whether to object returned should be of class \code{\link[sp]{SpatialPointsDataFrame}} for easier plotting with the \code{sp} package.  Default is \code{TRUE}.
#' @param dmethod The method used to decompose the covariance matrix for conditional simulation.  Valid options are \code{"chol"}, \code{"eigen"}, and \code{"svd"}.  The default is \code{"chol"}.
#' @param ... Currently unimplemented.
#' 
#' @return If \code{sp = TRUE}, then a \code{SpatialPointDataFrame} from the \code{sp} package is returned, with components including the prediction coordinates, the predicted responses \code{pred}, the mean square prediction error (\code{mspe}), the root mean square prediction error (\code{rmspe}), and the conditional realizations, if application (\code{sim.1}, \code{sim.2}, ...).  If \code{sp = FALSE}, then a list of class \code{gearKrige} is returned, with components \code{pred}, \code{mspe}, \code{rmspe}, and \code{sim}, if relevant.   
#' @author Joshua French
#' @importFrom stats predict delete.response rnorm model.matrix terms
#' @importFrom sp spDists
#' @examples 
#' 
#' # generate response
#' y = rnorm(10)
#' # generate coordinates
#' x1 = runif(10); x2 = runif(10)
#' 
#' # data frame for observed data
#' data = data.frame(y, x1, x2)
#' # newdata must have columns with prediction coordinates
#' newdata = data.frame(x1 = runif(5), x2 = runif(5))
#' 
#' # specify a standard covariance model
#' cmod = cmod.std(model = "exponential", psill = 1, 
#'                 r = 1)
#' 
#' # geolm for universal kriging
#' gearmod_uk = geolm(y ~ x1 + x2, data = data,
#'                  coordnames = c("x1", "x2"),
#'                  cmod = cmod)
#' # prediction for universal kriging, with conditional simulation
#' pred_uk = predict(gearmod_uk, newdata, nsim = 2)
#' 
#' # geolm for ordinary kriging
#' gearmod_ok = geolm(y ~ 1, data = data,
#'                  coordnames = c("x1", "x2"),
#'                  cmod = cmod)
#'# prediction for ordinary kriging
#' pred_ok = predict(gearmod_ok, newdata)
#' 
#' # geolm for simple kriging
#' gearmod_ok = geolm(y ~ 1, data = data,
#'                  coordnames = c("x1", "x2"),
#'                  cmod = cmod, mu = 1)
#'# prediction for simple kriging
#' pred_sk = predict(gearmod_ok, newdata)
#' 
#' @rdname predict.geolmStd
#' @export
predict.geolmStd = function(object, newdata, nsim = 0, 
                          sp = TRUE, dmethod = "chol", ...) 
{
  #check validity of arguments
  if(nsim < 0 || !is.finite(nsim))
  { stop("nsim should be a non-negative integer") }
  if(!is.element(dmethod, c("chol", "eigen", "svd")))
  {
    stop("Invalid dmethod type")
  }

  newcoords = as.matrix(newdata[,object$coordnames])
  vop = eval.cmod(object$cmod_evar0, 
                d = sp::spDists(object$coords, as.matrix(newcoords), 
                                longlat = object$longlat))
  
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
    vp = eval.cmod(object$cmod_evar0, 
                 d = sp::spDists(as.matrix(newcoords), 
                                 longlat = object$longlat))
    
    n = nrow(object$coords); m = nrow(newcoords)
    # correct mean if simple kriging
    mu = 0
    if(!is.null(object$mu)) mu = object$mu
    
    newsim = decomp.cov(rbind(
      cbind(object$v - diag(object$vediag),vop),
      cbind(t(vop), vp)), method=dmethod) %*% 
      matrix(rnorm((n+m)*nsim), ncol = nsim) + 
      mu
    
    # slightly different algorithm for simple kriging
    # update various object to include simulated data
    object$y = cbind(object$y, newsim[1:n,] + 
                matrix(stats::rnorm(n*nsim, 
                       sd = sqrt(object$vediag)),
                       nrow = n, ncol = nsim))
    
    #object$cholviy = backsolve(object$cholv, object$y, 
    #                            transpose = TRUE, upper.tri = TRUE)
    if(!is.null(object$mu))
    {
      object$cholviresid = 
        backsolve(object$cholv, object$y - mu, 
                  transpose = TRUE)
    }else
    {
      object$coeff = object$map2coeff %*% object$y
      object$cholviresid = backsolve(object$cholv, object$y - object$x %*% object$coeff, transpose = TRUE)
    }
  }
  
  cholvivop = forwardsolve(object$cholv, vop, 
                           transpose = TRUE, 
                           upper.tri = TRUE)
  
  # simple kriging
  if(!is.null(object$mu))
  {
    pred = object$mu + crossprod(cholvivop, object$cholviresid)
    
    # W = solve(object$v, vop) # weights matrix for prediction
    # list of predicted values and associated mse
    # pred2 = object$mu + crossprod(W, object$y - object$mu)
    # range(pred - pred2)
    mspe = object$cmod_evar0$psill + object$cmod_evar0$fvar - 
      colSums(cholvivop^2)
    mspe[mspe < 0] = 0 # numerical imprecision correction
    
    # mspe2 = colSums((object$v %*% W) * W) - 2 * colSums(W * vop) + object$cmod_evar0$psill + object$cmod_evar0$fvar
    # range(mspe - mspe2)
  }else
  {
    #compute kriging weights
    # 
    # vix <- solve(object$v, object$x)
    # xtvix <- crossprod(vix, object$x)
    # W <- solve(object$v, vop - object$x %*% solve(xtvix, crossprod(vix, vop) - t(newx)))
    # pred2 = crossprod(W, object$y)
    
    # cholviresid = forwardsolve(object$cholv, 
    #                            object$y - x %*% object$coeff,
    #                            transpose = TRUE, 
    #                            upper.tri = TRUE)
    pred = c(newx %*% object$coeff) + crossprod(cholvivop, object$cholviresid)
    
    # range(pred - pred2)
    # A = newx - crossprod(cholvivop, object$cholvix)
    # mspe2 = diag(vp) - diag(t(vop) %*% solve(object$v, vop)) +  diag(A %*% solve(xtvix, t(A)))

    mspe = object$cmod_evar0$psill + object$cmod_evar0$fvar - 
            colSums(cholvivop^2) + 
            colSums(forwardsolve(object$cholxtvix, 
                    t(newx - crossprod(cholvivop, 
                                       object$cholvix)), 
                    transpose = TRUE, upper.tri = TRUE)^2)
    mspe[mspe < 0] = 0 # numerical imprecision correction
    # range(mspe - mspe2)
    # list of objects to return
 }
#   # generate conditional realizations if nsim > 0
#   if(nsim > 0)
#   { 
#     out$sim = newsim[-(1:n),] + (pred[,1] - pred[,-1])
#   }
  
  if(sp)
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

