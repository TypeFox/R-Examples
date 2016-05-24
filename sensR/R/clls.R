clls <-
  function(location, scale, data, weights, start,..., subset,
           na.action, contrasts = NULL, Hess = FALSE, model = TRUE,
           method = c("logistic", "probit", "cloglog", "cauchit"))
{
    .Deprecated("clm", package="sensR",
                msg="'clls' is deprecated: use 'clm' from package 'ordinal' instead")
  logit <- function(p) log(p/(1 - p))
  fmin <- function(beta) {
    theta <- beta[pc + 1:q]
    ctheta <- c(-100, theta, 100)
    eta <- moffset
    if (pc > 0)
      eta <- eta + drop(x %*% beta[1:pc])
    esigma <- exp(Moffset)
    if(k > 0)
      esigma <- esigma * exp(drop(z %*% beta[pc+q + 1:k]))
    pr <- pfun((ctheta[y + 1] - eta)/esigma ) -
      pfun((ctheta[y] - eta)/esigma)
    if (all(pr > 0))
      -2*sum(wt * log(pr))
    else Inf
  }
  gmin <- function(beta) {
    theta <- beta[pc + 1:q]
    ctheta <- c(-100, theta, 100)
    eta <- moffset
    sigma <- exp(Moffset)
    if(k > 0)
      sigma <- sigma * exp(drop(z %*% beta[pc+q + 1:k]))
    if (pc > 0)
      eta <- eta + drop(x %*% beta[1:pc])
    m1 <- ctheta[y + 1] - eta
    m2 <- ctheta[y] - eta
    M1 <- m1/sigma
    M2 <- m2/sigma
    pr <- pfun(M1) - pfun(M2)
    p1 <- dfun(M1)
    p2 <- dfun(M2)
    g3 <- if(k > 0)
      t(z) %*% (wt * (m1*p1 - m2*p2)/(pr * sigma))
    else numeric(0)
    g1 <- if (pc > 0)
      t(x) %*% (wt * (p1 - p2)/(pr * sigma))
    else numeric(0)
    xx <- .polrY1 * p1 - .polrY2 * p2
    g2 <- -t(xx) %*% (wt/(pr * sigma))
    g2 <- t(g2)
    if (all(pr > 0))
      c(g1, g2, g3)
    else rep(NA, pc + q + k)
  }

  ### warning("This function is no longer supported.\nThe user is adviced to use clm() from package ordinal instead.")
  m <- match.call(expand.dots = FALSE)
  if(missing(location))
    stop("Model needs a specification of the location")
  if(missing(scale))
    m$scale <- ~1 ## No model for the scale is assumed
  method <- match.arg(method)
  pfun <- switch(method, logistic = plogis, probit = pnorm,
                 cloglog = pgumbel, cauchit = pcauchy)
  dfun <- switch(method, logistic = dlogis, probit = dnorm,
                 cloglog = dgumbel, cauchit = dcauchy)
  if (is.matrix(eval.parent(m$data)))
    m$data <- as.data.frame(data)
  m$start <- m$Hess <- m$model <- m$... <- m$method <- NULL
  m[[1]] <- as.name("model.frame")
  M <- m
  m$scale <- NULL
  M$location <- NULL
  names(m)[names(m) == "location"] <- "formula"
  names(M)[names(M) == "scale"] <- "formula"
  m <- eval.parent(m)
  Termsm <- attr(m, "terms")
  x <- model.matrix(Termsm, m, contrasts)
  xint <- match("(Intercept)", colnames(x), nomatch = 0)
  n <- nrow(x)
  pc <- ncol(x)
  wt <- model.weights(m)
  if(!length(wt)) {
    wt <- rep(1, n)
    M$weights <- wt ## Trick to avoid warning message in
    ##'eval.parent(M)' when scale part is omitted.
  }
  M <- eval.parent(M)
  TermsM <- attr(M, "terms")
  z <- model.matrix(TermsM, M, contrasts)
  zint <- match("(Intercept)", colnames(z), nomatch = 0)
  k <- ncol(z)
  cons <- attr(x, "contrasts")
  if (xint > 0) {
    x <- x[, -xint, drop = FALSE]
    pc <- pc - 1
  }
  else warning("an intercept is needed and assumed")
  if (zint > 0) {
    z <- z[, -zint, drop = FALSE]
    k <- k - 1
  }
  if(k > 0 && n != nrow(z))
    stop("Model needs same dataset in location and scale")
  moffset <- model.offset(m)
  if(length(moffset) <= 1)
    moffset <- rep(0, n)
  Moffset <- model.offset(M)
  if(length(Moffset) <= 1)
    Moffset <- rep(0, n)
  y <- model.response(m)
  if(!is.factor(y))
    stop("response must be a factor")
  lev <- levels(y)
  if (length(lev) <= 2)
    stop("response must have 3 or more levels")
  y <- unclass(y) # numeric levels
  q <- length(lev) - 1 # No. estimable intercepts
  Y <- matrix(0, n, q)
  .polrY1 <- col(Y) == y
  .polrY2 <- col(Y) == y - 1

  ## Get starting values:
  if(missing(start)) {
    ## try logistic/probit regression on 'middle' cut
    q1 <- length(lev) %/% 2
    y1 <- (y > q1)
    X <- cbind(Intercept = rep(1, n), x)
    fit <-
      switch(method,
             "logistic"= glm.fit(X, y1, wt, family = binomial(), offset = moffset),
             "probit" = glm.fit(X, y1, wt, family = binomial("probit"), offset = moffset),
             ## this is deliberate, a better starting point
             "cloglog" = glm.fit(X, y1, wt, family = binomial("probit"), offset = moffset),
             "cauchit" = glm.fit(X, y1, wt, family = binomial("cauchit"), offset = moffset))
    if(!fit$converged)
      stop("attempt to find suitable starting values failed")
    coefs <- fit$coefficients
    if(any(is.na(coefs))) {
      warning("design appears to be rank-deficient, so dropping some coefs")
      keep <- names(coefs)[!is.na(coefs)]
      coefs <- coefs[keep]
      x <- x[, keep[-1], drop = FALSE]
      pc <- ncol(x)
    }
    spacing <- logit((1:q)/(q+1)) # just a guess
    if(method != "logit") spacing <- spacing/1.7
    thetas <- -coefs[1] + spacing - spacing[q1]
    start <- c(coefs[-1], thetas, rep(0, ncol(z)))
  } else if(length(start) == pc + q + k) {
    if(k > 0)
      start[(pc+q+1):(pc+q+k)] <- log(start[(pc+q+1):(pc+q+k)])
  }
  else
    stop("'start' is not of the correct length")
  ## Optimize the log-likelihood function:
  fit <- optim(start, fn = fmin,  gr=gmin,
                method = "BFGS", hessian = Hess, ...)
  ## Extract parameters:
  beta <- fit$par[seq_len(pc)]
  theta <- fit$par[pc + 1:q]
  tau <- fit$par[pc + q + seq_len(k)]
  sigma <- exp(tau)
  deviance <- fit$value
  niter <- c(f.evals = fit$counts[1], g.evals = fit$counts[2])
  names(theta) <- paste(lev[-length(lev)], lev[-1], sep = "|")
  if(k > 0) {
    names(sigma) <- colnames(z)
    names(tau) <- colnames(z)
    Sigma <- exp(drop(z %*% tau))
  }
  else
    Sigma <- 1
  if (pc > 0) {
    names(beta) <- colnames(x)
    eta <- drop(x %*% beta)
  }
  else
    eta <- rep(0, n)
  cumpr <- matrix(pfun((matrix(theta, n, q, byrow = TRUE) - eta)/Sigma), ,
                  ncol = q)
  fitted <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
  fitted.case <- fitted[col(fitted) == y]

  dimnames(fitted) <- list(row.names(m), lev)
  res <- list(coefficients = c(beta, theta, sigma),
              beta = beta, theta = theta, sigma = sigma,
              tau = tau, deviance = deviance,
              fitted = fitted.case,
              fitted.values = fitted, lev = lev, terms.location =
              Termsm, terms.scale = TermsM,
              df.residual = sum(wt) - pc - q - k, edf = pc + q + k,
              n = sum(wt), nobs = sum(wt),
              call = match.call(),
              convergence = fit$convergence,
              niter = niter)
  if(Hess) {
    dn <- c(names(beta), names(theta), names(sigma))
    H <- fit$hessian
    dimnames(H) <- list(dn, dn)
    res$Hessian <- H
  }
  if(model) {
    res$location <- m
    res$scale <- M
  }
  res$na.action <- attr(m, "na.action")
  res$contrasts <- cons
  res$xlevels <- .getXlevels(Termsm, m)
  res$zlevels <- .getXlevels(TermsM, M)
  class(res) <- "clls"
  res
}

print.clls <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control=NULL)
    }
    if(length(x$beta)) {
        cat("\nLocation coefficients:\n")
        print(x$beta, ...)
    } else {
        cat("\nNo location coefficients\n")
    }
    if(length(x$sigma)) {
        cat("\nScale coefficients:\n")
        print(x$sigma, ...)
    } else {
        cat("\nNo Scale coefficients\n")
    }

    cat("\nIntercepts:\n")
    print(x$theta, ...)
    cat("\nResidual Deviance:", format(x$deviance, nsmall=2), "\n")
    cat("AIC:", format(x$deviance + 2*x$edf, nsmall=2), "\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
    if(x$convergence > 0)
        cat("Warning: did not converge as iteration limit reached\n")
    invisible(x)
}

vcov.clls <- function(object, ...)
{
  if(is.null(object$Hessian)) {
    message("\nRe-fitting to get Hessian\n")
    utils::flush.console()
    object <- update(object, Hess=TRUE, start = coef(object))
  }

  structure(solve(object$Hessian), dimnames =
            dimnames(object$Hessian))
}

summary.clls <- function(object, digits = max(3, .Options$digits - 3),
                         correlation = FALSE, ...)
{
    cc <- c(object$beta, object$theta, object$tau)
    pc <- length(object$beta)
    q <- length(object$theta)
    k <- length(object$tau)
    coef <- matrix(0, pc+q+k, 3, dimnames =
                   list(names(cc), c("Value", "Std. Error", "z value")))
    coef[, 1] <- cc
    vc <- vcov(object)
    coef[, 2] <- sd <- sqrt(diag(vc))
    coef[, 3] <- coef[, 1]/coef[, 2]
    object$coefficients <- coef
    object$pc <- pc
    object$q <- q
    object$k <- k
    object$digits <- digits

    if(correlation)
        object$correlation <- (vc/sd)/rep(sd, rep(pc+q+k, pc+q+k))
    class(object) <- "summary.clls"
    object
}

print.summary.clls <- function(x, digits = x$digits, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control=NULL)
    }
    coef <- format(round(x$coefficients, digits=digits))
    pc <- x$pc
    q <- x$q
    k <- x$k
    if(pc > 0) {
        cat("\nLocation coefficients:\n")
        print(coef[seq_len(pc), , drop=FALSE],
              quote = FALSE, ...)
    } else {
        cat("\nNo location coefficients\n")
    }
    if(k > 0) {
      cat("\nScale coefficients:\n")
      print(coef[(pc+q+1):(pc+q+k), , drop=FALSE],
            quote = FALSE, ...)
###       cc <- x$coefficients[(pc+q+1):(pc+q+k), , drop = FALSE]
###       cc <- cbind(exp(cc[, 1]), cc)
###       colnames(cc)[1:2] <- c("sigma", "log(sigma)")
###       cc <- format(round(cc, digits=digits))
###       print(cc, quote = FALSE, ...)
    } else {
      cat("\nNo scale coefficients\n")
    }
    cat("\nIntercepts:\n")
    print(coef[(pc+1):(pc+q), , drop=FALSE], quote = FALSE, ...)
    cat("\nResidual Deviance:", format(x$deviance, nsmall=2), "\n")
    cat("AIC:", format(x$deviance + 2*x$edf, nsmall=2), "\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
    if(!is.null(correl <- x$correlation)) {
        cat("\nCorrelation of Coefficients:\n")
        ll <- lower.tri(correl)
        correl[ll] <- format(round(correl[ll], digits))
        correl[!ll] <- ""
        print(correl[-1, -ncol(correl)], quote = FALSE, ...)
    }
    invisible(x)
}

anova.clls <- function (object, ..., test = c("Chisq", "none"))
{
  test <- match.arg(test)
  dots <- list(...)
  if (length(dots) == 0)
    stop('anova is not implemented for a single "clls" object')
  mlist <- list(object, ...)
  nt <- length(mlist)
  dflis <- sapply(mlist, function(x) x$df.residual)
  s <- order(dflis, decreasing = TRUE)
  mlist <- mlist[s]
  if (any(!sapply(mlist, inherits, "clls")))
    stop('not all objects are of class "clls"')
  ns <- sapply(mlist, function(x) length(x$fitted.values))
  if(any(ns != ns[1]))
    stop("models were not all fitted to the same size of dataset")
  rsp <- unique(sapply(mlist, function(x) {
                       tmp <- x$terms.location
                       class(tmp) <- "formula"
                       paste(tmp[2]) } ))
  mds <- sapply(mlist, function(x) {
                tmp1 <- x$terms.location
                tmp2 <- x$terms.scale
                class(tmp1) <- class(tmp2) <- "formula"
                paste(tmp1[3], "|", tmp2[2]) } )
  dfs <- dflis[s]
  lls <- sapply(mlist, function(x) deviance(x))
  tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
  df <- c(NA, -diff(dfs))
  x2 <- c(NA, -diff(lls))
  pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
  out <- data.frame(Model = mds, Resid.df = dfs, Deviance = lls,
                    Test = tss, Df = df, LRtest = x2, Prob = pr)
  names(out) <- c("Model", "Resid. df", "Resid. Dev", "Test",
                  "   Df", "LR stat.", "Pr(Chi)")
  if (test == "none") out <- out[, 1:6]
  class(out) <- c("Anova", "data.frame")
  attr(out, "heading") <-
    c("Likelihood ratio tests of cumulative link location-scale models\n",
      paste("Response:", rsp))
  out
}

logLik.clls <- function(object, ...)
  structure(-0.5 * object$deviance, df = object$edf, class = "logLik")

extractAIC.clls <- function(fit, scale = 0, k = 2, ...)
{
    edf <- fit$edf
    c(edf, deviance(fit) + k * edf)
}

pgumbel <- function(q, loc = 0, scale = 1, lower.tail = TRUE)
{
    q <- (q - loc)/scale
    p <- exp(-exp(-q))
    if (!lower.tail) 1 - p else p
}

dgumbel <- function (x, loc = 0, scale = 1, log = FALSE)
{
    d <- log(1/scale) - x - exp(-x)
    if (!log) exp(d) else d
}

