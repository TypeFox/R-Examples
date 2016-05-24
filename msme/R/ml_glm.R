ml_glm <- function(formula,
                   data,
                   family,
                   link,
                   offset = 0,
                   start = NULL,
                   verbose = FALSE,
                   ...) {

### Handle the input
  mf <- model.frame(formula, data)
  y <- model.response(mf, "numeric")

### Prepare model infrastructure
  class(y) <- c(family, link, "expFamily")
  X <- model.matrix(formula, data = data)

### Check for missing data.  Stop if any.
  if (any(is.na(cbind(y, X)))) stop("Some data are missing!")

### Initialize the search, if needed
  if (is.null(start))  start <- kickStart(y, X, offset)
  
### Maximize the joint log likelihood
  fit <- maximize(start, Sjll, X, y, offset, ...)

### Check for optim failure and report and stop
  if (verbose | fit$convergence > 0)  print(fit)

### Extract and compute quantities of interest
  beta.hat <- fit$par
  se.beta.hat <- sqrt(diag(solve(-fit$hessian)))
  residuals <- devianceResiduals(y, beta.hat, X, offset, ...)

### Fit null model and determine null deviance 
  fit.null <- maximize(mean(y), Sjll, 1, y, offset, ...)
  null.deviance <-
    sum(devianceResiduals(y, fit.null$par, 1, offset, ...)^2)

### Report the results, with the needs of print.glm in mind
  results <- list(fit = fit,
                  X = X,
                  y = y,
                  call = match.call(),
                  obs = length(y),
                  df.null = length(y) - 1,
                  df.residual = length(y) - length(beta.hat),
                  deviance = sum(residuals^2),    
                  null.deviance = null.deviance, 
                  residuals = residuals,
                  coefficients = beta.hat,
                  se.beta.hat = se.beta.hat,
                  aic = - 2 * fit$val + 2 * length(beta.hat),
                  i = fit$counts[1])

### Use (new) msme class and glm class
  class(results) <- c("msme","glm")
  return(results)
}
