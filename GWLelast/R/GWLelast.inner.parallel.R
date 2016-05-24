#' Inner part of fitting GWLelast with parallel cores
#'
#' @param x Covariates.
#' @param y Outcome binary variable.
#' @param coords 2 columns matrix including "longitude" and "latitude".
#' @param W Weight matrix.
#' @param alpha The elasticnet mixing parameter [0,1] in glmnet package.
#' @param lambda Optional user-supplied lambda sequence in glmnet package.
#' @param nlambda The number of lambda values in glmnet package.
#' @return model Fitted model at location i.
#' @return cv.error Cross validation error.

GWLelast.inner.parallel = function(x=x, y = y, coords = coords, W = W, alpha = 1, lambda = lambda, nlambda = nlambda) {
  
  model = list()
  cv.error = list()
  GWLelast.object = list()
  
  i = 1
  
  result = foreach(i = 1:dim(x)[1],.packages=c('glmnet'), .errorhandling='remove') %dopar% {        
    w = W[i, -i]
    
    xi = as.matrix(x[-i,])
    yi = as.matrix(y[-i])
    
    model = glmnet(x = xi, y = yi, weights = w, family = "binomial", alpha = alpha, nlambda = nlambda, lambda = lambda)
    predictions = predict(model, newx = t(x[i,]), type = "class") 
    cv.error = abs(as.numeric(predictions) - y[i])
    
    return(list(model=model, cv.error=cv.error))
  }
  
  GWLelast.object[["model"]] = sapply(result,function(x){list(x$model)})
  GWLelast.object[["cv.error"]] = sapply(result,function(x){list(x$cv.error)})
  
  return(GWLelast.object)
}