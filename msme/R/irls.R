### IRLS --- ### IRLS ---  GLM ESTIMATION FUNCTION  :   -----irls.r-----
### Hilbe, J., and A. Robinson, Methods of Statistical Model Estimation
###    Chapman & Hall/CRC (2011), Chapter 4     10 July, 2011 revision

irls <- function(formula, data, family, link,
                 tol = 1e-6,
                 offset = 0,
                 m = 1,
                 a = 1,
                 verbose = 0) {

### Prepare the model components as previously
  mf <- model.frame(formula, data)
  y <- model.response(mf, "numeric")
  X <- model.matrix(formula, data = data)

### Check for missing values 
  if (any(is.na(cbind(y, X)))) stop("Some data are missing.")

### Arrange the class information
  class(y) <- c(family, link, "expFamily")

### Establish a start point
  mu <- initialize(y, m)
  eta <- linkFn(mu, m, a)
  dev <- devIRLS(y, mu, m, a)
  deltad <- 1
  i <- 0

### IRLS loop (class is preserved by these arithmetic ops)
  while (abs(deltad) > tol ) {
    w <-  1 / (variance(mu, m, a) * lPrime(mu, m, a)^2)
    z <- eta + (y - mu) * lPrime(mu, m, a) - offset
    mod <- lm(z ~ X - 1, weights = w)
    eta <- mod$fit + offset
    mu <- unlink(y, eta, m, a)
    dev.old <- dev
    dev <- devIRLS(y, mu, m, a)
    deltad <- dev - dev.old
    i <- i + 1
    if(verbose > 0) cat(i, coef(mod), deltad, "\n")
  }

### Post-estimation statistics
  df.residual <- summary(mod)$df[2]
  pearson.chi2 <- sum((y - mu)^2 /
                  variance(mu, m, a)) / df.residual
  ll <- sum(jllm(y, mu, m, a))

### Return a rich object --- allows use of print.glm
  result <-
    list(coefficients = coef(mod),
         se.beta.hat = sqrt(diag(summary(mod)$cov.unscaled)),
         model = mod,
         call = match.call(),
         nobs = length(y),
         eta = eta,
         mu = mu,
         df.residual = df.residual,
         df.null = length(y) - 1,
         deviance = dev,
         null.deviance = NA,
         p.dispersion = pearson.chi2,
         pearson = pearson.chi2 * df.residual,
         loglik = ll,
         family = list(family = family),
         X = X,
         i = i,
         residuals = devianceResids(y, mu, m, a),
         aic = -2 * ll + 2 * summary(mod)$df[1])
  class(result) <- c("msme","glm")
  return(result)
}
