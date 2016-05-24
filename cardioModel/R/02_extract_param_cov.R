extract.param.cov <- function(model, parameters) {
  ret <- data.frame()
  # Extract the coefficients or fixed effects using the right
  # function
  if (any(c("lme", "nlme") %in% class(model))) {
    tmp.coef <- fixef(model)
  } else {
    tmp.coef <- coef(model)
  }
  tmp.vcov <- vcov(model)
  # Iterate over all parameters for values and covariances
  for (i in 1:length(parameters)) {
    # Extract the parameters
    ret[1,paste("value", parameters[i], sep=".")] <-
      tmp.coef[parameters[i]]
    # Iterate from the current to the last parameter for the covariances
    for (j in i:length(parameters)) {
      ret[1,paste("cov", parameters[i], parameters[j], sep=".")] <-
        tmp.vcov[parameters[i], parameters[j]]
    }
  }
  ret
}
