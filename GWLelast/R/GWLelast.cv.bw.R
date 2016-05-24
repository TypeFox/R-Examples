#' Cross validation for geographically weighted logistic elastic net regression
#'
#' @param x Covariates.
#' @param y Outcome binary variable.
#' @param coords 2 columns matrix including "longitude" and "latitude".
#' @param alpha The elasticnet mixing parameter [0,1] in glmnet package.
#' @param lambda Optional user-supplied lambda sequence in glmnet package.
#' @param nlambda The number of lambda values in glmnet package.
#' @param gweight geographical kernel function in spgwr package.
#' @param longlat Indicate if the coords parameter are sperically calculated.
#' @param bw bandwidth of geographical kernel function.
#' @param D Distance matrix.
#' @param Parallel Calculate the model with multi core or not.
#' @return cv.error Cross validation error.

GWLelast.cv.bw = function(x = x, y = y, coords = coords, alpha = 1, lambda = lambda, nlambda = nlambda, gweight = gweight, longlat = longlat, bw = bw, D = D, Parallel = Parallel) {
  
  cat(paste("Initial bandwidth: ", bw, "\n", sep = ""))
  
  GWLelast.model = GWLelast.est(x = x, y = y, coords= coords, alpha = alpha, lambda = lambda, 
                                      nlambda = nlambda, gweight = gweight, 
                                      longlat = longlat, bw = bw, D = D, Parallel = Parallel)
  
  cv.error = sum(sapply(GWLelast.model[["cv.error"]], sum))
  
  cat(paste("CV error: ", cv.error, "\n", sep = ""))
  
  return(cv.error)
}
