lqa.default <- function (x, y, family = gaussian (), penalty = NULL, method = "lqa.update2", weights = rep (1, nobs), start = NULL, etastart = NULL, mustart = NULL, offset = rep (0, nobs), control = lqa.control (), intercept = TRUE, standardize = TRUE, ...)
{
    call <- match.call ()

### Check for exponential family and link function:
### -----------------------------------------------

    if (is.character (family)) 
      family <- get (family, mode = "function", envir = parent.frame ())

    if (is.function (family)) 
      family <- family ()

    if (is.null (family$family)) 
    {
      print (family)
      stop ("'family' not recognized")
    }


### Check for quadratic penalty:
### ---------------------------- 

  if (! (method == "nng.update"))
  {

    if (is.null (penalty))
      stop ("penalty not specified \n")

    if (is.character (penalty)) 
      penalty <- get (penalty, mode = "function", envir = parent.frame ())

    if (is.function (penalty)) 
      penalty <- penalty ()

    if (is.null (penalty$penalty))
    {
      print (penalty)
      stop ("'penalty' not recognized")
    }
  }


### Check for existence of 'method':
### --------------------------------

    if (is.null (method))
      stop ("method not specified")


### Check for column of ones in x if intercept == TRUE:
### ---------------------------------------------------

    if (intercept & (var (x[,1]) > control$var.eps))
      x <- cbind (1, x)


### Standardization:
### ----------------

    x <- as.matrix (x)
    xnames <- dimnames (x)[[2L]]

    ynames <- if (is.matrix (y)) 
                rownames(y)
              else 
                names(y)

    nobs <- nrow (x)   
    nvars <- ncol (x)    # number of coefficients in the predictor (including an intercept, if present)
    ones <- rep (1, nobs)
    mean.x <- drop (ones %*% x) / nobs       # computes the vector of means

    if (intercept)    # if an intercept is included in the model its corresponding element of mean.x is set to zero
      mean.x[1] <- 0  # (such that x[,1] is not getting centered (and also not standardized later on ...))

    x.std <- scale (x, mean.x, FALSE)   # centers the regressor matrix

    norm.x <- if (standardize)
              {
                 norm.x <- sqrt (drop (ones %*% (x.std^2)))   # computes the euclidean norm of the regressors
                 nosignal <- apply (x, 2, var) < control$var.eps
                 if (any (nosignal))    # modify norm.x for variables with too small a variance (e.g. the intercept)
                   norm.x[nosignal] <- 1

                 norm.x
              }
              else
                rep (1, nvars)

    x.std <- scale (x.std, FALSE, norm.x)  # standardizes the centered regressor matrix


### Call and get the (estimation) method: 
### -------------------------------------

    fit <- do.call (method, list (x = x.std, y = y, family = family, penalty = penalty, intercept = intercept, control = control, ...))


### Back-Transformation of estimated coefficients:
### ----------------------------------------------

    coef <- fit$coefficients
    
    if (intercept)
    {
       coef[1] <- coef[1] - sum (mean.x[-1] * coef[-1] / norm.x[-1])
       coef[-1] <- coef[-1] / norm.x[-1]
    }
    else
       coef <- coef / norm.x


### Computation of some important statistics:
### -----------------------------------------

# Remark: The predictors are identical no matter whether we use standardized or unstandardized values, hence all statistics
#         based on the predictor eta are also equal

    eta <- drop (x %*% coef)
    mu <- family$linkinv (eta)
    mu.eta.val <- family$mu.eta (eta)
    wt <- sqrt ((weights * mu.eta.val^2) / family$variance (mu))
    dev <- sum (family$dev.resids (y, mu, weights))
    wtdmu <- sum (weights * y) / sum (weights)
    nulldev <- sum (family$dev.resids (y, wtdmu, weights))
    n.ok <- nobs - sum (weights == 0)
    nulldf <- n.ok - as.integer (intercept)
    residuals <- (y - mu) / mu.eta.val

    xnames <- colnames (x)
    ynames <- names (y)
    names (residuals) <- names (mu) <- names (eta) <- names (weights) <- names (wt) <- ynames
    names (coef) <- xnames

    Amat <- fit$Amat
    Amat <- t (norm.x * t (norm.x * Amat))   # must be (quadratically) transformed in order to cope with the transformed parameter space

    if (is.null (fit$tr.H))
      stop ("quadpen.fit: Element 'tr.H' has not been returned from 'method'")

    model.aic <- dev + 2 * fit$tr.H
    model.bic <- dev + log (nobs) * fit$tr.H
    resdf <- n.ok - fit$tr.H

    dispersion <- ifelse (!((family$family == "binomial") | (family$family == "poisson")), sum ((y - mu)^2 / family$variance (mu)) / (nobs - fit$tr.H), 1)
    
    fit <- list (coefficients = coef, residuals = residuals, fitted.values = mu, family = family, penalty = penalty, 
linear.predictors = eta, deviance = dev, aic = model.aic, bic = model.bic, null.deviance = nulldev, n.iter = fit$stop.at, 
best.iter = fit$m.stop, weights = wt, prior.weights = weights, df.null = nulldf, df.residual = resdf, converged = fit$converged, mean.x = mean.x, 
norm.x = norm.x, Amat = Amat, method = method, rank = fit$tr.H, x = x, y = y, fit.obj = fit, call = call, dispersion = dispersion)
   
    class (fit) <- c ("lqa", "glm", "lm")
    fit
}
