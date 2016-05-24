betaCIContinuousExactTrans <-
function(object, level)
 {
  mod <- object$model
  coefs <- summary(mod)$coefficients
  beta <- coefs[2, 1]
  sdbeta <- coefs[2, 2]
  prob <- (1 + level) / 2
  if (any(family(mod)$family == c("binomial", "poisson")))
   z <- qnorm(prob) else  z <- qt(prob, df = mod$df.residual)
  betaCI <- beta + c(0, -1, 1) * z * sdbeta
  betaCI <- matrix(betaCI, nrow = 1, ncol = 3, byrow = TRUE)
  colnames(betaCI) <- c("beta", paste(50 * (1 + c(-1, 1) * level), "%"))
  betaCI
 }
