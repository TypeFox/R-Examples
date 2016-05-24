ml_g <- function(formula, data) {

### Prepare the data, relying on the formula class for 
### handling model specification
  mf <- model.frame(formula, data)
  y <- model.response(mf, "numeric")
  X <- model.matrix(formula, data = data)

### Check for missing data.  Stop if any.
  if (any(is.na(cbind(y, X)))) stop("Some data are missing.")

### Declare the joint log likelihood function 
  jll.normal <- function(params, X, y) {
     p <- length(params)
     beta <- params[-p]
     sigma <- exp(params[p])
     linpred <- X %*% beta  
     sum(dnorm(y, mean = linpred, sd = sigma, log = TRUE))
     }

### Initialize the search
  ls.reg <- lm(y ~ X - 1)
  beta.hat.ls <- coef(ls.reg)
  sigma.hat.ls <- sd(residuals(ls.reg))
  start <- c(beta.hat.ls, sigma.hat.ls)

### Maximize the joint log likelihood
  fit <- optim(start,           
               jll.normal,
               X = X,
               y = y,
               control = list(
                 fnscale = -1,
                 maxit = 10000),
               hessian = TRUE
               )

### Check for optim failure and report and stop
  if (fit$convergence > 0) {
    print(fit)
    stop("optim failed to converge!")
  }

### Post-processing 
  beta.hat <- fit$par
  se.beta.hat <- sqrt(diag(solve(-fit$hessian)))

### Reporting
  results <- list(fit = fit,
                  X = X,
                  y = y,
                  call = match.call(),
                  beta.hat = beta.hat,
                  se.beta.hat = se.beta.hat,
                  sigma.hat = exp(beta.hat[length(beta.hat)]))

### Prepare for S3 deployment (see next Section!)
  class(results) <- c("ml_g_fit","lm")
  return(results)
}
