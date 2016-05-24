ml_glm3 <- function(formula,
                    data,
                    family,
                    link,
                    offset = 0,
                    start = NULL,
                    verbose = FALSE,
                    ...) {
  mf <- model.frame(formula, data)
  y <- model.response(mf, "numeric")
  class(y) <- c(family, link, "expFamily")
  X <- model.matrix(formula, data = data)
  if (any(is.na(cbind(y, X)))) stop("Some data are missing!")
  if (is.null(start))  start <- kickStart(y, X, offset)
  fit <- maximize(start, Sjll, X, y, offset, ...)
  if (verbose | fit$convergence > 0)  print(fit)
  beta.hat <- fit$par
  se.beta.hat <- sqrt(diag(solve(-fit$hessian)))
  residuals <- devianceResiduals(y, beta.hat, X, offset, ...)
  results <- list(fit = fit,
                  X = X,
                  y = y,
                  call = match.call(),
                  obs = length(y),
                  df.null = length(y) - 1,
                  df.residual = length(y) - length(beta.hat),
                  deviance = sum(residuals^2),    
                  null.deviance = NA,
                  residuals = residuals,
                  coefficients = beta.hat,
                  se.beta.hat = se.beta.hat,
                  aic = - 2 * fit$val + 2 * length(beta.hat),
                  i = fit$counts[1])
  class(results) <- c("msme","glm")
  return(results)
}
