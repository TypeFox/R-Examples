# NB1 maximum likelihood function J Hilbe and A Robinson 11Apr 2010, 10Jul 2011

ml.nb1 <- function(formula, data, offset = 0, start = NULL, verbose = FALSE) {
### Prepare the data, relying on the formula class for 
### handling model specification
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  nb1X <- model.matrix(formula, data = data)
### Check for missing data.  Stop if any.
  if (any(is.na(cbind(y, nb1X)))) stop("Some data are missing.")
### Declare the joint log likelihood function 
  jll.nb1 <- function(b.hat, X, y) {
    a.hat <- b.hat[1]
    xb.hat <- X %*% b.hat[-1] + offset
    mu.hat <- exp(xb.hat)
    r.hat <- (1/a.hat) * mu.hat    
    sum(dnbinom(y,
                size = r.hat,
                mu = mu.hat,
                log = TRUE))
    }
### Initialize the search
  if (is.null(start))
    start <- c(0.5, coef(lm.fit(nb1X, log(y+1))))
### Maximize the joint log likelihood
  fit <- optim(start,           
               jll.nb1,
               X = nb1X,
               y = y,
               control = list(
                 fnscale = -1,
                 maxit = 10000),
               hessian = TRUE
               )
### Check for optim failure and report and stop
  if (verbose | fit$convergence > 0) print(fit)
### Compute model for null deviance
  fit.null <- optim(start[1],           
                    jll.nb1,
                    X = 1,
                    y = y,
                    method = "BFGS",
                    control = list(
                      fnscale = -1)
                    )
  beta.hat <- fit$par
  se.beta.hat <- sqrt(diag(solve(-fit$hessian)))
  y.hat <- exp(nb1X %*% beta.hat + offset)
  residuals <- NA
  results <- list(fit = fit,
                  X = nb1X,
                  y = y,
                  y.hat = y.hat,
                  call = match.call(),
                  obs = length(y),
                  df.null = length(y) - 1,
                  df.residual = length(y) - length(start),
                  deviance = - 2 * fit$val,
                  null.deviance = - 2 * fit.null$val,
                  residuals = residuals,
                  loglike = fit$val, 
                  aic = - 2 * fit$val + 2 * length(start),
                  coefficients = beta.hat,
                  se.beta.hat = se.beta.hat)
  class(results) <- c("msme","glm")
  return(results)
}





