# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

fastLasso <- function(x, y, lambda, subset = NULL, normalize = TRUE, 
                      intercept = TRUE, eps = .Machine$double.eps, 
                      use.Gram = TRUE, drop = TRUE, raw = FALSE) {
  # initializations
  intercept <- isTRUE(intercept)
  use.Gram <- isTRUE(use.Gram)
  drop <- isTRUE(drop)
  raw <- isTRUE(raw)
  # compute lasso
  if(raw) {
    # call C++ function
    fit <- callBackend("R_testLasso", R_x=x, R_y=y, R_lambda=lambda, 
                       R_initial=seq_along(y), R_intercept=intercept, 
                       R_eps=eps, R_useGram=use.Gram)
    
    # prepare object for raw lasso fit
    coef <- fit$coefficients
    res <- fit$residuals
    if(drop) {
      # drop the dimension of the components
      coef <- drop(coef)
      res <- drop(res)
    }
    center <- mean(res)
    scale <- sqrt(mean((res-center)^2))
    fit <- list(best=fit$indices, coefficients=coef, residuals=res, 
                objective=fit$crit, center=center, scale=scale)
  } else {
    # check subset
    if(is.null(subset)) {
      useSubset <- FALSE
      subset <- integer()
    } else useSubset <- TRUE
    # call C++ function
    fit <- callBackend("R_fastLasso", R_x=x, R_y=y, R_lambda=lambda, 
                       R_useSubset=useSubset, R_subset=subset, 
                       R_normalize=normalize, R_intercept=intercept, 
                       R_eps=eps, R_useGram=use.Gram)
    if(drop) fit <- lapply(fit, drop)  # drop the dimension of the components
  }
  # return lasso fit
  fit
}
