discrimR <-
  function(formula, data, weights, cluster, start, subset, na.action,
           contrasts = NULL, hess = FALSE, ranef = FALSE, zi = FALSE,
           method = c("duotrio", "probit", "threeAFC", "triangle",
             "twoAFC"), ...)
{
  ## ...: arguments to optim, such as control=list(trace=TRUE)
  ## start: optional starting values

### Function to minimized:
  fminNM <- function(beta) {
    for(j in 1:nl) { 
      func1 <- function(u) {
        x <- X[ind==j, ,drop=FALSE] ## within cluster
        y <- Y[ind==j]
        w <- wt[ind==j]
        eta <- offset[ind==j] + x %*% beta[1:nb] + u ## linear pred.
        p <- link.inv(eta)
        b <- exp(beta[nb + ns]) ## exp(sigma)
        d <- dnorm(u/b)/b ## normal density
        fp0[j] <<- p0^(w * y) * (1 - p0)^(w * (1 - y)) * pnorm(-beta[1]/b)
        prod( p^(w * y) * (1 - p)^(w * (1 - y)) ) * d
      }
      h[j] <<- integrate(Vectorize(func1), -beta[1], Inf)$value
      f[j] <<- fp0[j] + h[j]
    }
    if (all(f > 0)) 
      -2*sum(log(f)) # - minus twice negative log likelihood
    else Inf
  }
  m <- match.call(expand.dots = FALSE)
  m$start <- m$hess <- m$method <- m$... <- m$ranef <- m$zi <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  X <- model.matrix(Terms, m, contrasts)
  n <- nrow(X)
  cons <- attr(X, "contrasts")
  wt <- model.weights(m)
  if (!length(wt)) 
    wt <- rep(1, n)
  offset <- model.offset(m)
  if (length(offset) <= 1) 
    offset <- rep(0, n)
  Y <- model.response(m)
  if (NCOL(Y) == 2) {
    n <- Y[, 1] + Y[, 2]
    Y <- ifelse(n == 0, 0, Y[, 1]/n)
    wt <- wt * n
  }
  stopifnot(all(wt >= 0))
  method <- match.arg(method)
  link.inv <- switch(method,
                     duotrio = duotrio()$linkinv,
                     probit = pnorm,
                     threeAFC = threeAFC()$linkinv,
                     triangle = triangle()$linkinv,
                     twoAFC = twoAFC()$linkinv,
                     )

  ## Parameters to fminNM:
  nb <- NCOL(X); nb
  ns <- 1; ns ## 1
  nl <- nlevels(m$"(cluster)"); nl
  ind <- c(unclass(m$"(cluster)")); ind ## indicater for clusters
  f <- h <- g <- fp0 <- double(nl)
  p0 <- ifelse(method %in% c("triangle", "threeAFC"), 1/3, .5)
  
  ## Starting values:
  suppressWarnings( {
    if(missing(start)) {
      gf <- glm.fit(X, Y, family=binomial("probit")) # added "" to probit 1.4-6
      b1 <- as.vector(gf$coef); b1
      start <- c(b1, 1)
    }
  })
  if(length(start) != nb + ns)
    stop("'start' is not of correct length")
  if(start[nb + ns] <= 0)
    stop("starting value for sigma needs to possitive")
  start <- c(start[1:nb], log(start[nb + ns]))
  
  ## Optimization
  fit <- optim(par=start, fn=fminNM, method="BFGS", 
               hessian=hess, ...) 

  ## Output:
  fpar <- fit$par[seq_len(nb)]; fpar
  rpar <- c(fit$par[nb+ns], exp(fit$par[nb+ns])); rpar
  se <- if(hess) sqrt(diag(solve(fit$hessian)))
  else NULL
  deviance <- fit$value; deviance
  p <- 1 - pnorm(-fit$par[1]/exp(fit$par[nb + ns]))
  resid.p <- resid <- fitted <- NULL
  
  ## Ranef:
  if(ranef) {
    for(j in 1:nl) {
      func2 <- function(u) {
        x <- X[ind==j, ,drop=FALSE] ## within cluster
        y <- Y[ind==j]
        w <- wt[ind==j]
        eta <- offset[ind==j] + x %*% fit$par[1:nb] + u ## linear pred.
        p <- link.inv(eta)
        b <- exp(fit$par[nb + ns]) ## exp(sigma)
        d <- dnorm(u/b)/b ## normal density
        prod( p^(w * y) * (1 - p)^(w * (1 - y)) ) * (fit$par[1] + u) * d
      }
      g[j] <- integrate(Vectorize(func2), -fit$par[1], Inf)$value
    }
    ranef <- g/f
    fitted <- link.inv(ranef)
    resid <- sign(Y - fitted) * sqrt(-2*log(f))
    V <- fitted * (n - fitted)/n
    resid.p <- (Y - fitted)/sqrt(V)
  }
  else
    {
    ranef <- NULL
    g <- NULL
  }
  if(zi) {
    zi <- h / f
  }
  else {
    zi <- NULL
  }
      
  list(fpar=fpar, rpar=rpar, deviance=deviance, se=se, #par = fit$par, 
       convergence=fit$convergence, lli = f, ranef = ranef, zi
       = zi, p = p, fitted = fitted, Y = Y, call = match.call())
}
