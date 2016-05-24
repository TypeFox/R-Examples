ocv <- function(mod){       # leave-one-out cross-validation
  y <- mod$model[,1] 
  b <- coef(mod)
  A <- hatvalues(mod)
  ocv <- mean((y-fitted(mod))^2/(1-A)^2)
  return(ocv)
}
