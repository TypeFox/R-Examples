################################
#### Compositional data ridge regression tuning of the
#### lambda and the alpha parameters via k-fold cross validation
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
#### References:Tsagris M. (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean journal of statistics 6(2): 47-57
################################

alfa.ridge <- function(y, x, a, lambda, B = 1, xnew = NULL) {
  ## y is dependent univariate variable. It can be matrix or vector
  ## x are compositional data, the covariates
  ## a is the value of a for the alpha-transformation
  ## lambda is the ridge regularization parameter
  ## if lambda=0, the classical multivariate regression is implemented
  ## B is for bootstrap estimation of the standard errors of the betas
  ## if pred is TRUE it means that you want to predict new y
  ## but if xnew is x (by default), the pred is not important
  ## the pred is important if xnew is not x
  z <- alfa(x, a, h = TRUE)$aff ## apply the alpha-transformation
  mod <- ridge.reg(y, z, lambda, B = B, xnew = xnew)
  mod ## ridge regression with the alpha-transformed compositional data
}
