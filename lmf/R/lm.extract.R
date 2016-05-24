lm.extract <-
function(formula, data, na.action = na.exclude)
{
  #"mod" is a fitted linear model
  mod <- lm(formula = formula, data = data, na.action = na.action)
  #Alphas
  ret <- list(ajt = mod$coefficients)
  #Residuals
  ret$res <- residuals(mod)
  #Fitted values
  ret$fit <- predict(mod)
  #Degrees of freedom
  ret$dof <- df.residual(mod)
  #Residual standard error
  ret$sigma.djt <- sqrt(deviance(mod)/ret$dof)
  #Unscaled covariance matrices (the variance-covariance
  #matrix is given by "cov.unscaled * sigma^2")
  ret$Ajt.us <- vcov(mod)/ret$sigma.djt^2
  # Add if sentence here to account for cases of no variation in y (Wj) Ajt.us should be all zeros and not NaN then.
  if(ret$sigma.djt == 0){
    ret$Ajt.us <- vcov(mod)
  }
  #Leverage
  ret$leverage <- lm.influence(mod, do.coef = FALSE)$hat
  #Cooks distance
  ret$cook <- cooks.distance(mod, sd = ret$sigma.djt, res = ret$res)
  #Output
  ret
}
