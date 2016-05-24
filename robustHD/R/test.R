# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## for testing purposes

# test object-based lasso
testLasso <- function(x, y, lambda, initial, intercept = TRUE, 
                      eps = .Machine$double.eps, use.Gram = TRUE) {
  callBackend("R_testLasso", R_x=x, R_y=y, R_lambda=lambda, R_initial=initial, 
              R_intercept=intercept, R_eps=eps, R_useGram=use.Gram)
}

# test object-based C-step
testCStep <- function(x, y, lambda, subset, h, intercept = TRUE, 
                      tol = .Machine$double.eps^0.5, eps = .Machine$double.eps, 
                      use.Gram = TRUE) {
  callBackend("R_testCStep", R_x=x, R_y=y, R_lambda=lambda, R_subset=subset, 
              R_intercept=intercept, R_tol=tol, R_eps=eps, R_useGram=use.Gram)
}

# test keeping best subsets
testKeepBest <- function(subsetMat, crits, nkeep) {
  callBackend("R_testKeepBest", R_subsetMat=subsetMat, 
              R_crits=crits, R_nkeep=nkeep)
}
