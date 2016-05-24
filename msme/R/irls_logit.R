irls_logit <- function(formula, data, tol = 0.000001) {# arguments

## Set up the model components
    mf <- model.frame(formula, data)        # define model frame as mf
    y <- model.response(mf, "numeric")      # set model response as y
    X <- model.matrix(formula, data = data) # predictors in matrix X

## Check for missing values; stop if any.
    if (any(is.na(cbind(y, X)))) stop("Some data are missing.")

## Initialize mu, eta, the deviance, etc.
    mu <- rep(mean(y), length(y))
    eta <- log(mu/(1-mu))                             
    dev <- 2 * sum(y*log(1/mu) +
           (1 - y) * log(1/(1-mu)))
    deltad <- 1                                       
    i <- 1                                            

## Loop through the IRLS algorithm
    while (abs(deltad) > tol ) {                  # IRLS loop begin
        w <-  mu * (1-mu)                         # weight
        z <- eta + (y - mu)/w                     # working response
        mod <- lm(z ~ X-1, weights = w)           # weighted regression
        eta <- mod$fit                            # linear predictor
        mu <- 1/(1+exp(-eta))                     # fitted value
        dev.old <- dev                            
        dev <- 2 * sum(y * log(1/mu) +
                       (1 - y) * log(1/(1 - mu))) # deviance
        deltad <- dev - dev.old                   # change
        cat(i, coef(mod), deltad, "\n")           # iteration log
        i <- i + 1                                # iterate
    }

## Build some post-estimation statistics
    df.residual <- summary(mod)$df[2]
    pearson.chi2 <- sum((y - mu)^2 / (mu * (1 - mu))) / df.residual

## Return a compact result
    result <- list(coefficients = coef(mod),
                   se.beta.hat = sqrt(diag(summary(mod)$cov.unscaled)))
  return(result)
}
