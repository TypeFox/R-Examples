tsglm.check <- function(fit){
  if(!("tsglm" %in% class(fit))) stop("Object does not have class 'tsglm'")
  listnames <- c("coefficients", "start", "residuals", "fitted.values", "linear.predictors", "response", "logLik", "score", "info.matrix", "info.matrix_corrected", "call", "n_obs", "n_eff", "ts", "model", "xreg", "distr", "distrcoefs", "sigmasq")
  if(!(all(listnames %in% names(fit)))) stop("Object does not have the required named list elements.")
  #include more checks wether it is a proper object of class 'tsglm'
}