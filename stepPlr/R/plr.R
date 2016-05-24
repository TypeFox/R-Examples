plr <- function(x, y, weights = rep(1,length(y)),
                offset.subset = NULL, offset.coefficients = NULL,
                lambda = 1e-4, cp = 'bic')
  {
    this.call <- match.call()
    nobs <- length(y)
    mu0 <- sum(y * weights) / sum(weights)
    if (all(x == 1)) {
      if (is.vector(x))
        x <- matrix(1, nobs, 1)
      p <- 1
      level <- NULL
      eta <- b <- log(mu0 / (1 - mu0))
      mu <- rep(mu0, nobs)
      eta <- rep(eta, nobs)
    } else {
      if (!is.numeric(x)) {
        imat <- get.imat(x)
        x <- imat$imat
        level <- imat$level
      } else {
        level <- NULL
      }
      x <- cbind(1, x)
      p <- ncol(x)
      if (!is.null(offset.subset)) {
        offset.subset <- offset.subset + 1
        offset.v <- x[ , offset.subset, drop = FALSE] %*% offset.coefficients
        has.offset <- 1
        z <- c(as.vector(x[ , -offset.subset]), y, weights,
               as.vector(offset.v), has.offset)
      } else {
        has.offset <- 0
        z <- c(as.vector(x), y, weights, has.offset)
      }
      q <- p - length(offset.subset)
      sol <- .Fortran('solveplr',
                      x = as.double(rep(0, q)),
                      n = as.integer(q),
                      z = as.double(z),
                      lenz = as.integer(length(z)),
                      nobs = as.integer(nobs),
                      lam = as.double(lambda),
                      status = as.integer(-1))
      if (sol$status != 0)
        cat('\nConvergence warning in plr:', sol$status, '\n')
      if (has.offset == 1) {
        b <- rep(NA, p)
        b[offset.subset] <- offset.coefficients
        b[-offset.subset] <- sol$x
      } else {
        b <- sol$x
      }
      eta <- drop(x %*% b)
      mu <- 1 / (1 + exp(-eta))
    }
    w <- mu * (1 - mu)
    lam <- c(0, rep(lambda, (p - 1)))
    ddf0 <- t(x) %*% (x * weights * w)
    ddf <- ddf0 + 2 * diag(p) * lam
    ddf.i <- solve(ddf)
    cov <- ddf.i %*% ddf0
    df <- sum(diag(cov))
    cov <- cov %*% ddf.i
    devresids <- binomial()$dev.resids
    dev <- sum(devresids(y, mu, weights))
    null.dev <- sum(devresids(y, mu0, weights))
    if (cp == 'bic') {
      cp <- log(nobs)
    } else if (cp == 'aic') {
      cp <- 2
    }
    score <- dev + cp * df
    if (p == 1) {
      xnames <- 'Intercept'
    } else if (is.null(dimnames(x)[[2]])) {
      xnames <- c('Intercept', paste('x', seq(p - 1), sep = ''))
    } else {
      xnames <- c('Intercept', dimnames(x)[[2]][2:p])
    }
    names(b) <- xnames
    dimnames(cov) <- list(xnames, xnames)
    names(mu) <- names(eta) <- seq(nobs)
    object <- list(coefficients = b, covariance = cov, deviance = dev,
                   null.deviance = null.dev, df = df, score = score,
                   nobs = nobs, cp = cp, fitted.values = mu,
                   linear.predictors = eta, level = level, call = this.call)
    class(object) <- 'plr'
    return(object)
  }

predict.plr <- function(object, newx = NULL,
                        type = c('link', 'response', 'class'), ...)
  {
    type <- match.arg(type)
    if (is.null(newx)) {
      if (type == 'link') {
        pred <- object$linear.predictors
      } else if (type == 'response') {
        pred <- object$fitted.values
      } else {
        pred <- ifelse(object$fitted.values > 0.5, 1, 0)
      }
    }
    else {
      if (!is.null(object$level))
        newx <- get.imat(newx, object$level)$imat
      if (!all(newx==1)) {
        x <- cbind(1, newx)
      } else {
        x <- newx
      }
      if (is.vector(x))
        x <- matrix(x)
      pred <- drop(x %*% object$coefficients)
      if (type == 'response') {
        pred <- 1 / (1 + exp(-pred))
      } else if (type=='class') {
        pred <- ifelse(pred > 0, 1, 0)
      }
      names(pred) <- seq(length(pred))
    }
    return(pred)
  }

summary.plr <- function(object, ...)
  {
    cat('\nCall:\n')
    print(object$call)
    cat('\nCoefficients:\n')
    coef <- round(object$coefficients, digits=5)
    stderr <- round(sqrt(diag(object$covariance)), digits=5)
    zvalue <- round(coef / stderr, digits=3)
    pvalue <- round(2*(pnorm(-abs(zvalue))), digits=3)
    xnames <- names(coef)
    junk <- cbind(coef, stderr, zvalue, pvalue)
    dimnames(junk) <- list(names(coef),
                           c('Estimate','Std.Error','z value','Pr(>|z|)'))
    print(junk)
    cat('\n    Null deviance:', round(object$null.deviance, digits = 2),
        'on', object$nobs - 1,'degrees of freedom\n')
    cat('Residual deviance:', round(object$deviance, digits = 2),
        'on', round(object$nobs - object$df, digits = 2),
        'degrees of freedom\n')
    cat('            Score: deviance +', round(object$cp, digits=1),
        '* df =', round(object$score, digits = 2), '\n')
    return(invisible(junk))
  }

print.plr <- function(x, ...)
  {
    cat('\nCall:\n')
    print(x$call)
    cat('\nCoefficients:\n')
    print(round(x$coefficients, digits=5))
    cat('\n    Null deviance:', round(x$null.deviance, digits = 2),
        'on', x$nobs - 1,'degrees of freedom\n')
    cat('Residual deviance:', round(x$deviance, digits = 2),
        'on', round(x$nobs - x$df, digits = 2), 'degrees of freedom\n')
    cat('            Score: deviance +', round(x$cp, digits = 1),
        '* df =', round(x$score, digits = 2),'\n')
  }
