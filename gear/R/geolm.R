#' Linear model for geostatistical data.
#' 
#' \code{geolm} creates a geostatistical linear model object of the appropriate class based on the arguments, especially the \code{cmod} arguments.
#' 
#' \code{formula} should be specified after the form \code{y ~ x1 + x2}, where \code{y} is the response variable and \code{x1} and \code{x2} are the covariates of interest.  If \code{mu} is provided, the variables to the right of \code{~} are ignored.
#' 
#' @param formula An object of class formula.  See Details.
#' @param data A data frame containing the response, covariates, and location coordinates.
#' @param coordnames A vector of length 2 with the names of the columns in \code{data} containing the coordinates, e.g., \code{c("long", "lat")}.
#' @param cmod A covariance model object obtained from one of the \code{cmod.*} functions, e.g., \code{\link[gear]{cmod.std}}.
#' @param vmod A semivariance model object obtained from one of the \code{vmod.*} functions.  Not currently implemented.
#' @param weights An optional vector of weights for the errors to be used in the fitting process.  A vector that is proportional to the reciprocal variances of the errors, i.e., errors are assumed to be uncorrelated with variances \code{evar/weights}.  Default is \code{NULL}, meaning that the weights are uniformly 1.
#' @param longlat A logical value.  Default is \code{FALSE}.  If \code{TRUE}, Great Circle distances (WGS84 ellipsoid) are calcualted between locations.  Otherwise, Euclidean. 
#' @param mu A single numeric value indicating the consant mean of the spatial process if simple kriging is desired.  Default is \code{NULL}, meaning that ordinary or universal kriging should be used.
#' 
#' @return Returns a \code{geolm} object.
#' 
#' @author Joshua French
#' @importFrom sp spDists
#' @importFrom stats model.frame model.response model.matrix
#' @export
#' @examples 
#' data = data.frame(y = rnorm(10), x1 = runif(10), 
#'                  x2 = runif(10))
#' cmod = cmod.std(model = "exponential", psill = 1, 
#'                  r = 1)
#' gearmod = geolm(y ~ x1, data = data,
#'                  coordnames = c("x1", "x2"),
#'                  cmod = cmod)

geolm = function(formula, data, coordnames, 
                  cmod = NULL, vmod = NULL, 
                  weights = NULL, longlat = FALSE,
                  mu = NULL)
{
  check.args.geolm(formula, data, coordnames, cmod,
                    vmod, weights, longlat, mu)

  # create model frame
  mf = stats::model.frame(formula, data)
  # extract response
  y = c(stats::model.response(mf))
  n = length(y) # number of observations
  if(is.null(mu)) # if not doing simple kriging
  {
    # create design matrix, remove all attributes
    # to make it a basic matrix
    x = stats::model.matrix(formula, data = mf)
    dimnames(x) = NULL
  }else
  {
    x = NULL
  }

  # create default weights, if NULL
  if(is.null(weights)) weights = rep(1, n)
    
  if(!is.null(cmod))
  {
    if(is.element(class(cmod), c("cmodStd", "cmodMan")))
    {
      vediag = cmod$evar/weights
      # create covariance matrix for observed data
      coords = as.matrix(data[,coordnames])
      
      if(class(cmod) == "cmodStd")
      {
        cmod_evar0 = cmod
        cmod_evar0$evar = 0
        
        v = eval.cmod(cmod_evar0, 
                    d = sp::spDists(coords, longlat = longlat)) + 
                    diag(vediag)
        
      }else # it's a cmodMan
      {
        cmod_evar0 = NULL
        v = cmod$v
      }
      
      cholv = chol(v)
      
      if(!is.null(mu))
      {
        vix = NULL
        xtvix = NULL
        resid = y - mu
        coeff = NULL
        cholvix = NULL
        cholxtvix = cholxtvix = NULL 
        vcov = NULL
        map2coeff = map2coeff = NULL
      }else
      {
        ###compute matrix products for future use
        # vix <- solve(v, x)
        cholvix = backsolve(cholv, x, transpose = TRUE)
        vix = forwardsolve(cholv, cholvix, upper.tri = TRUE)
        # range(vix - vix2)
        # xtvix2 <- crossprod(vix, x)
        xtvix = crossprod(cholvix)
        cholxtvix = chol(xtvix)
        vcov = chol2inv(cholxtvix)
        # range(xtvix - xtvix2)
        #compute gls estimates of regression coefficients
        map2coeff = forwardsolve(cholxtvix, 
                                 backsolve(cholxtvix, t(vix), 
                                           transpose = TRUE), 
                                 upper.tri = TRUE)
        
        coeff = map2coeff %*% y
        # coeff2 <- solve(xtvix, crossprod(vix, y))
        # range(coeff - coeff2)
        resid = y - x %*% coeff
      }
      
      cholviresid = backsolve(cholv, resid, transpose = TRUE)
      # viresid = solve(v, resid)
      # viresid2 = forwardsolve(cholv, cholviresid, upper.tri = TRUE)
      # range(viresid - viresid2)
      # compute log likelihood
      loglik = -n/2*log(2*pi) - 
        #determinant(v, log = TRUE)$mod/2 - 
        sum(2 * log(diag(cholv)))/2 - 
        crossprod(cholviresid)[1,1]/2
        # crossprod(resid, solve(v, resid))[1,1]/2
      
      which_class = "geolmStd"
      if(class(cmod) == "cmodMan") which_class = "geolmMan"
      
      return(structure(list(y = y, x = x, loglik = loglik,
                 resid = resid, coeff = coeff, mu = mu, 
                 v = v, 
                 cholv = cholv,
                 cholvix = cholvix, # cholviy = cholviy, 
                 # vix = vix,
                 cholviresid = cholviresid,
                 cholxtvix = cholxtvix, 
                 vcov = vcov,
                 # xtvix = xtvix,
                 map2coeff = map2coeff, 
                 formula = formula, coordnames = coordnames,
                 cmod_evar0 = cmod_evar0,
                 weights = weights, 
                 vediag = vediag,
                 evar = cmod$evar,
                 coords = coords,
                 longlat = longlat), class = which_class))
    }
  }
}

check.args.geolm = function(formula, data, coordnames, 
                             cmod, vmod, 
                             weights, longlat, mu)
{

  if(class(formula) != "formula") 
    stop("formula must be a formula object or NULL")
  if(length(formula) != 3) stop("formula should have a single response to the left of the '~' and the covariates to the right")
  if(length(formula[[2]]) != 1) stop("There should only be a single response.")
  if(!is.data.frame(data)) stop("data should be a data frame")
  df_names = names(data)
  if(min(is.element(coordnames, df_names)) == 0) stop("The names in coordnames are not found in the column names of the data dataframe.")
  if(is.null(cmod) & is.null(vmod))
  {
    stop("Either cmod or vmod must be specified")
  }
  if(!is.null(cmod) & !is.null(vmod))
  {
    stop("Only one of cmod or vmod can be specified")
  }
  if(!is.null(cmod))
  {
    if(!is.element(class(cmod), c("cmodMan", "cmodStd")))
    {
      stop("cmod must be of class cmodMan or cmodStd")
    }
  }
  if(!is.null(weights))
  {
    if(length(weights) != nrow(data))
    {
      stop("length(weights) != nrow(data)")
    }
  }
  if(!is.logical(longlat))
  {
    stop("longlat must be a logical value")
  }
  if(!is.null(mu))
  {
    if(length(mu) != 1) stop("mu must be a vector of length 1")
    if(!is.numeric(mu)) stop("mu must be a numeric value")
  }
}


