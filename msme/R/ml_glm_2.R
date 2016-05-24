ml_glm2 <- function(formula1,
                    formula2 = ~1,
                    data,
                    family,
                    mean.link,
                    scale.link,
                    offset = 0,
                    start = NULL,
                    verbose = FALSE) {

### Handle the input
  mf <- model.frame(formula1, data)
  y <- model.response(mf, "numeric")

### Prepare model infrastructure
  class(y) <- c(family, mean.link, scale.link, "expFamily")
  X1 <- model.matrix(formula1, data = data)
  X2 <- model.matrix(formula2, data = data)
  colnames(X2) <- paste(colnames(X2), "_s", sep="")
  p <- ncol(X1)
  X <- cbind(X1, X2)

### Check for missing data.  Stop if any.
  if (any(is.na(cbind(y, X)))) stop("Some data are missing!")

### Initialize the search 
  if (is.null(start)) {
    start <- c(kickStart(y, X1, offset),
               1, # Shameless hack
               rep(0, ncol(X) - p - 1))
    names(start) <- c(colnames(X1), colnames(X2))
  }

### Maximize the joint log likelihood
  fit <- maximize(start, Sjll2, X, y, offset, p = p)

### Check for optim failure and report and stop
  if (verbose | fit$convergence > 0)  print(fit)

### Extract and compute quantities of interest
  beta.hat <- fit$par
  se.beta.hat <- sqrt(diag(solve(-fit$hessian)))
  residuals <- devianceResiduals2(y, beta.hat, X, p, offset)

#### Deviance residuals for null 
  fit.null <- maximize(c(mean(y), 1), 
                        Sjll2,
                        X[,c(1,p+1)],  y, offset, p = 1)
  null.deviance <-
    sum(devianceResiduals2(y,
                          c(fit.null$par[1], fit$par[p+1]),
                          X[, c(1,p+1)],
                          1,
                          offset)^2)

### Report the results, with the needs of print.glm in mind
  results <- list(fit = fit,
                  loglike = fit$val,
                  X = X,
                  y = y,
                  p = p,
                  call = match.call(),
                  obs = length(y),
                  df.null = length(y) - 2,
                  df.residual = length(y) - length(beta.hat),
                  deviance = sum(residuals^2),    
                  null.deviance = null.deviance, 
                  residuals = residuals,
                  coefficients = beta.hat,
                  se.beta.hat = se.beta.hat,
                  aic = - 2 * fit$val + 2 * length(beta.hat),
                  offset = offset,
                  i = fit$counts[1])
  class(results) <- c("msme","glm")
  return(results)
}
