#' Inner part of fitting GWLelast without parallel cores
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

GWLelast.inner = function(x=x, y = y, coords = coords, W = W, lambda = lambda, alpha = 1, nlambda = nlambda) {
  
  model = list()
  cv.error = list()
  GWLelast.inner = list()
  
  for(i in 1:dim(x)[1]){        
    w = W[i, -i]
    
    xi = as.matrix(x[-i,])
    yi = as.matrix(y[-i])
    
    model[[i]] = glmnet(x = xi, y = yi, weights = w, family = "binomial", alpha = alpha, nlambda = nlambda, lambda = lambda)
    predictions = predict(model, newx = t(x[i,]), type = "class") 
    cv.error[[i]] = abs(as.numeric(predictions) - y[i])
    
  }
  
  GWLelast.inner[["model"]] = model
  GWLelast.inner[["cv.error"]] = cv.error
  
  return(GWLelast.inner)
}