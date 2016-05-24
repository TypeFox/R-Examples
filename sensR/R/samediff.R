samediff <-
  function (nsamesame, ndiffsame, nsamediff, ndiffdiff, VCOV = TRUE)
{
  call <- m <- match.call()
  m[[1]] <- as.name("list")
  m <- eval.parent(m) # evaluate the *list* of arguments
  ss <- m[[1]]; ds <- m[[2]]; sd <- m[[3]]; dd <- m[[4]]
  data <- unlist(m)
  for(i in data) {
    if(i != trunc(i) | i < 0)
      stop("data have to be non-negative integers")
  }
  if(length(data[data == 0]) > 2)
    stop("Not enough information in data,",
         " when more than two entries are zero")
  
  vcov <- array(NA, dim = c(2,2))
  name <- c("tau", "delta")
  dimnames(vcov) <- list(name, name)
  se <- rep(NA, 2)
  conv <- NULL
  
  ## If the fraction of same answers to same-samples is larger than the
  ## fraction of same-answers to diff-samples, simple estimation of
  ## delta may proceed: 
  isdelta <- ifelse(ss/(ds + ss) > sd/(sd + dd), TRUE, FALSE)

  ## Estimate tau, delta and their covariance matrix:
  if(isdelta && all(data != 0)) {
    tau <- Tau(ss, ds)
    psd <- sd/(sd + dd)
    delta <- uniroot(f = Delta, lower = 0, upper = 10, tau = tau,
                      psd = psd, tol = 1e-6)$root
    i11 <- I11(tau, delta, ss, ds, sd, dd)
    i12 <- I12(tau, delta, sd, dd)
    i22 <- I22(tau, delta, sd, dd)
    vcov <- solve(cbind(c(i11, i12), c(i12, i22)))
    dimnames(vcov) <- list(name, name)
    se <- sqrt(diag(vcov))
    case <- 0
    logLik <- llSameDiff(tau, delta, ss, ds, sd, dd)
  }
  else if(ss == 0 && sd == 0) {
###    warning("No information on delta") ## parse this warning?
    tau <- 0
    delta <- NA
    logLik <- 0
    case <- 0.1
  }
  else if(ds == 0 && dd == 0) {
###    warning("No information on delta") ## parse this warning?
    tau <- Inf
    delta <- NA
    case <- 1
    logLik <- 0
  }
  else if(ds == 0 && sd == 0) {
    tau <- Inf
    delta <- Inf
    case <- 1.2
    logLik <- 0
  }
  else if(ss == 0 && ds == 0){
    tau <- Inf
    delta <- NA
    case <- 1.12
    tmp <- optim(c(5, 5), fn = llds0, sd = sd, dd = dd,
                 method = "L-BFGS-B", hessian = FALSE,
                 lower = c(1e-4, 1e-4), upper = c(100, 100), 
                 control = list(fnscale = -1))
    logLik <- tmp$value
    conv <- tmp$convergence
  }
  else if(sd == 0 && dd == 0) {
    delta <- NA
    tmp <- optim(1, fn = llDeltaInf, ss = ss, ds = ds,
                 method = "BFGS", hessian = VCOV,
                 control = list(fnscale = -1))
    tau <- tmp$par
    if(VCOV) {
      vcov[1, 1] <- -(tmp$hessian)^(-1)
      se <- c(sqrt(vcov[1, 1]), NA)
    }
    conv <- tmp$convergence
    case <- 1.3
    logLik <- llDeltaInf(tau, ss, ds)
  }
  else if(ds == 0) {
    tau <- Inf
    delta <- Inf
    case <- 1.22
    tmp <- optim(c(5, 5), fn = llds0, sd = sd, dd = dd,
                 method = "L-BFGS-B", hessian = FALSE,
                 lower = c(1e-4, 1e-4), upper = c(100, 100), 
                 control = list(fnscale = -1))
    logLik <- tmp$value
    conv <- tmp$convergence
  }
  else if(sd == 0) {
    delta <- Inf
    tmp <- optim(1, fn = llDeltaInf, ss = ss, ds = ds,
                 method = "BFGS", hessian = VCOV,
                 control = list(fnscale = -1))
    tau <- tmp$par
    if(VCOV) {
      vcov[1, 1] <- -(tmp$hessian)^(-1)
      se <- c(sqrt(vcov[1, 1]), NA)
    }
    conv <- tmp$convergence
    case <- 3
    logLik <- llDeltaInf(tau, ss, ds)
  }
  else if(dd == 0 || ss == 0 || !isdelta) {
    delta <- 0
    tmp <- optim(1, fn = llSameDiff, delta = 1e-4, ss = ss, ds = ds,
                 sd = sd, dd = dd, method = "BFGS",
                 hessian = VCOV,
                 control = list(fnscale = -1))
    tau <- tmp$par
    if(VCOV) {
      vcov[1, 1] <- -(tmp$hessian)^(-1)
      se <- c(sqrt(vcov[1, 1]), NA)
    }
    conv <- tmp$convergence
    case <- 2
    logLik <- llSameDiff(tau, 1e-4, ss, ds, sd, dd)
  }
  else
    stop("This should never occur! Please contact package maintainer!") 
  
  coef <- c(tau, delta)
  names(se) <- names(coef) <- name  
  fit <- list(coef = coef, se = se, vcov = vcov, data = data,
              test = "same-different",
              call = call, convergence = conv, case = case,
              logLik = logLik)
  class(fit) <- "samediff"
  fit
}

Tau <- function(ss, ds) { # compute tau
  frac <- (2 * ss + ds)/(2 * (ds + ss))
  tau <- sqrt(2) * qnorm(frac)
  tau
}

Delta <- function(d, tau, psd)  # compute d-prime
  pnorm((tau - d)/sqrt(2)) - pnorm(-(tau + d)/sqrt(2)) - psd


I11 <- function(tau, d, ss, ds, sd, dd) { ## 11-element of Fisher I-matrix
  E1 <- (dd + sd)^3/(2 * dd * sd)
  E2 <- (dnorm((tau - d)/sqrt(2)) + dnorm((-tau - d)/sqrt(2)))^2
  E3 <- 2*(ds + ss)^3/(ds * ss)
  E4 <- (dnorm(tau/sqrt(2)))^2
  E1 * E2 + E3 * E4
}

I11.tau <- function(tau, ss, ds) { ## Fisher information for tau
  E3 <- 2 * (ds + ss)^3/(ds * ss)
  E4 <- (dnorm(tau/sqrt(2)))^2
  E3 * E4
}

I12 <- function(tau, d, sd, dd) { ## 12 = 21-element of Fisher I-matrix
  E1 <- -(dd + sd)^3/(2 * dd * sd)
  E2 <- (dnorm((tau - d)/sqrt(2)))^2 - (dnorm((-tau - d)/sqrt(2)))^2
  E1 * E2
}

I22 <- function(tau, d, sd, dd) { ## 22-element of Fisher I-matrix
  E1 <- (dd + sd)^3/(2 * dd * sd)
  E2 <- (dnorm((tau - d)/sqrt(2)) - dnorm((-tau - d)/sqrt(2)))^2
  E1 * E2
}

llSameDiff <- function(tau, delta, ss, ds, sd, dd, verbose = FALSE) {
  if(any(c(tau, delta) <= 0)) return(-Inf)
  d <- delta
  sqrt.2 <- sqrt(2)
  Pss <- 2 * pnorm(tau/sqrt.2) - 1
  Pds <- 1 - Pss
  Psd <- pnorm((tau - d)/sqrt.2) - pnorm((-tau - d)/sqrt.2)
  Pdd <- 1 - Psd
  da <- c(ss, ds, sd, dd)
  pDa <- c(Pss, Pds, Psd, Pdd)
  da[pDa == 0] <- 0
  pDa[pDa == 0] <- 1
  ll <- ss * log(Pss) + ds * log(Pds) + sd * log(Psd) +
    dd * log(Pdd)
  if(verbose)
    cat(ll, tau, delta, Pss, Pds, Psd, Pdd, "\n")
  ll
}

llds0 <- function(par, sd, dd, verbose = FALSE) {
  if(any(par <= 0)) return(-Inf)
  d <- par[1]
  tau <- par[2]
  sqrt.2 <- sqrt(2)
  Psd <- pnorm((tau - d)/sqrt.2) - pnorm((-tau - d)/sqrt.2)
  Pdd <- 1 - Psd
  ll <- sd * log(Psd) + dd * log(Pdd)
  if(verbose)
    cat(ll, tau, d, Psd, Pdd, "\n")
  ll
}

llDeltaInf <- function(tau, ss, ds, verbose = FALSE) {
  if(any(tau <= 0)) return(-Inf)
  sqrt.2 <- sqrt(2)
  Pss <- 2 * pnorm(tau/sqrt.2) - 1
  Pds <- 1 - Pss
  ll <- ss * log(Pss) + ds * log(Pds)
  if(verbose)
    cat(ll, tau, Pss, Pds, "\n")
  ll
}

tauRootDeltaInf <- function(tau, ss, ds, pll.max, lim) {
### Solve this function for tau to obtain a confidence limit, when
### hat.delta == Inf
  pll <- llDeltaInf(tau, ss, ds)
  exp(pll - pll.max) - lim
}

ll2SameDiff <- function(delta, tau, ss, ds, sd, dd)
  llSameDiff(tau, delta, ss, ds, sd, dd) ## reverse arguments

print.samediff <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

vcov.samediff <- function(object, ...) {
  object$vcov
}

coef.samediff <- function(object, ...) {
  object$coef
}

logLik.samediff <- function(object, ...) {
  co <- coef(object)
  nobs <- sum(object$data)
  val <- object$logLik
  names(val) <- NULL
  attr(val, "nobs") <- nobs
  attr(val, "df") <- ifelse(is.na(co[2]), 1, 2)
### df = 2 even if parameters are known to be either 0 or Inf?
  class(val) <- "logLik"
  val
}

profile.samediff <-
  function(fitted, which = 1:2, max = 2, numpts = 100,
           max.delta = 10, max.tau = 10, ...)
{
  co <- fitted$coef
  da <- fitted$data
  
  ## Profiling tau:
  if(fitted$case == 1.12) {
    M <- max.tau
    ttau <- seq(1e-4, M, length = numpts) ## no evaluation at zero
    ll <- double(numpts)
    for(i in 1:numpts)
      ## optimize likelihood for delta at each of the tau-values.
      ll[i] <- pll.tau(delta = c(1e-2, M+2), tau = ttau[i], ss = da[1],
                       ds = da[2], sd = da[3], dd = da[4])
  }
  else if(is.na(co[2]) || co[1] == Inf && co[2] == Inf) {
    M <- max.tau
    ttau <- seq(1e-4, M, length = numpts) ## no evaluation at zero
    if(da[2] == 0 && da[3] == 0 || da[2]==0 && da[4]==0) 
      ll <- llDeltaInf(ttau, ss = da[1], ds = da[2])
    else {
      M <- max.tau
      ttau <- seq(1e-4, M, length = numpts) ## no evaluation at zero
      ll <- double(numpts)
      for(i in 1:numpts)
        ## optimize likelihood for delta at each of the tau-values.
        ll[i] <- pll.tau(delta = c(1e-2, M+2), tau = ttau[i], ss = da[1],
                         ds = da[2], sd = da[3], dd = da[4])
    }
  }
  else {
    if(co[1] == Inf)
      M <- max.tau
    else
      M <- co[1] + max
    ttau <- seq(1e-4, M, length = numpts) ## no evaluation at zero
    ll <- double(numpts)
    for(i in 1:numpts)
      ## optimize likelihood for delta at each of the tau-values.
      ll[i] <- pll.tau(delta = c(1e-2, M), tau = ttau[i], ss = da[1],
                       ds = da[2], sd = da[3], dd = da[4])
  }
  Data <- data.frame("tau" = ttau, "plTau" = ll)
  
  ## Profiling delta:
  if(!is.na(co[2])) {
    if(co[2] == Inf)
      M <- max.delta
    else
      M <- co[2] + max
    
    delta <- seq(1e-4, M, length = numpts)
    for(i in 1:numpts) {
      ll[i] <- pll.delta(tau = c(1e-4, M), delta = delta[i],
                         ss = da[1], ds =  da[2], sd = da[3], dd = da[4])
    }
    Data$delta <- delta
    Data$plDelta <- ll
  }
  attr(Data, "logLik") <- fitted$logLik
  class(Data) <- c("profile.samediff", "data.frame")
  ## Return:
  Data
}
  
plot.profile.samediff <-
  function(x, which = 1:nc, level = c(0.99, 0.95),
           fig = TRUE, ...)
{
  max.ll <- attr(x, "logLik")
  nc <- ncol(x)/2
  if(!all(which %in% 1:2) && length(which) <= 2)
    stop("Only 1 and 2 allowed in the vector 'which'",
         " with maximum two elements") 
  oldpar <- par(mfrow = c(1, length(which)))
  on.exit(par(oldpar))
  lim <- sapply(level, function(x) exp(-qchisq(x, df=1)/2) )
  x$nplTau <- exp(x[,2] - max.ll)
  x$nplDelta <- exp(x[,4] - max.ll)
  
  if(fig == TRUE) {
    if(any(which == 1)) {
      plot(x[,1], x$nplTau, type = "l", las = 1, ylim = c(0, 1), 
           xlab = expression(tau), 
           ylab = "Normalized Profile Likelihood",
         main = "")
      abline(h = lim, col = "grey")
    }
    if(any(which == 2)) {
      plot(x[,3], x$nplDelta, type = "l", las = 1, ylim = c(0, 1), 
           xlab = expression(delta),
           ylab = "Normalized Profile Likelihood",
           main = "")
      abline(h = lim, col = "grey")
    }
  }
  class(x) <- c("nProfile.samediff", "data.frame")
  invisible(x)
}

summary.samediff <- function(object, profile = TRUE, ...) {
  co <- coef(object)
  if(profile){
    CI <- confint(object, ...)
    pval <- pval.prof(object, ...)
  }
  else {
    CI <- CI.vcov(object, ...)
    pval <- pval.vcov(object, ...)
  }
  object$table <-
    cbind("Estimate" = co, "Std. Error" = SE.samediff(object), CI,
          "P-value" = pval) ## Names?
  object$AIC <- AIC(object)
  class(object) <- "summary.samediff"
  object
}

print.summary.samediff <-
  function(x, ..., digits = max(3, getOption("digits") - 3),
           signif.stars = getOption("show.signif.stars"))
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Coefficients\n")
  printCoefmat(x$table, digits = digits,
               signif.stars=signif.stars, tst.ind=integer(0),
               cs.ind=1:4,
               P.values=TRUE, has.Pvalue=TRUE,)
  cat("\nLog Likelihood:", round(x$logLik, digits),
      "\tAIC:", round(x$AIC, digits), "\n") 
  invisible(x)
}

pval.prof <- function(object, ...) {
  co <- coef(object)
  da <- object$data
  mll <- c(logLik(object))
  
  if(co[1] != 0) {
    if(co[1] == Inf && object$case %in% c(1, 1.2)) 
      PllTau <- llDeltaInf(1e-4, ss=da[1], ds=da[2])
    else 
      PllTau <- pll.tau(, 1e-4, ss=da[1], ds=da[2], sd=da[3],
                        dd=da[4])
    pval.tau <- 1 - pnorm(sqrt(2 * (mll - PllTau)))
  }
  else
    pval.tau <- 1

  if(is.na(co[2]))
    pval.delta <- NA
  else if(co[2] == 0)
    pval.delta <- 1
  else  {
    ## If delta > 0 and potentially Inf:
    PllDelta <-
      pll.delta(, 1e-4, ss=da[1], ds=da[2], sd=da[3], dd=da[4])
    pval.delta <- 1 - pnorm(sqrt(2 * (mll - PllDelta)))
  }
    
  pval <- c(pval.tau, pval.delta)
  pval
}

pval.vcov <- function(object, ...) {
  co <- coef(object)
  SE <- SE.samediff(object)
  wald.statistic <- co/SE
  pval <- 1 - pnorm(wald.statistic)
  pval
}

SE.samediff <- function(object, ...) {
  SE <- rep(NA, 2)
  var <- vcov(object)
  if(!is.na(var[1, 1]))
    SE[1] <- sqrt(var[1, 1])
  if(!is.na(var[2, 2]))
    SE[2] <- sqrt(var[2, 2])
  SE
}

CI.vcov <- function(object, level = .95, ...) {
  co <- coef(object)
  SE <- SE.samediff(object)
  lim <- qnorm((1 - level)/2, lower.tail = FALSE)
  CI <- matrix(rep(NA, 4), 2)
  dimnames(CI) <- list(c("", ""), c("Lower", "Upper"))

  if(!is.na(SE[1])) {
    CI[1, 1] <- co[1] - lim * SE[1]
    CI[1, 2] <- co[1] + lim * SE[1]
  }
  if(!is.na(SE[2])) {
    CI[2, 1] <- co[2] - lim * SE[2]
    CI[2, 2] <- co[2] + lim * SE[2]
  }

  CI
}

deltaRoot <- function(delta, ss, ds, sd, dd, pll.max, lim) {
### Solve this function for delta to obtain a confidence limit
  pll <- pll.delta( , delta, ss, ds, sd, dd)
  exp(pll - pll.max) - lim
}

tauRoot <- function(tau, ss, ds, sd, dd, pll.max, lim) {
### Solve this function for tau to obtain a confidence limit
  pll <- pll.tau( , tau, ss, ds, sd, dd)
  exp(pll - pll.max) - lim
}

pll.delta <- function(tau = c(1e-4, 10), delta, ss, ds, sd, dd) {
### Profile log-likelihood for delta
  optimize(llSameDiff, tau, delta = delta, ss = ss,
           ds = ds, sd = sd, dd = dd,
           maximum = TRUE)$objective
}

pll.tau <- function(delta = c(1e-4, 10), tau, ss, ds, sd, dd) {
### Profile log-likelihood for tau
  optimize(ll2SameDiff, delta, tau = tau, ss = ss,
           ds = ds, sd = sd, dd = dd, maximum = TRUE)$objective
}

confint.samediff <-
  function(object, parm = c("tau", "delta"), level = 0.95,
           max = c(10, 10), ...)
{
  co <- coef(object)
  da <- object$data
  ss <- da[1]; ds <- da[2]; sd <- da[3]; dd <- da[4]
  lim <- exp(-qchisq(level, df = 1)/2)
  lim2 <- exp(-qchisq(level, df = .5)/2)
  CI <- matrix(rep(c(0, Inf), each = 2), 2)
  lev.names <- paste(level, "% ", c("tau", "delta"), sep = "")
  dimnames(CI) <- list(lev.names, c("Lower", "Upper"))
  co[co == 0] <- 1e-4
  if(is.na(co[2])){ 
    parm <- ifelse("tau" %in% parm, "tau", "")
    CI[2, ] <- rep(NA, 2)
  }
  pll.max <- c(object$logLik)
  
  if("tau" %in% parm) {
    ## TAU:  
    if(object$case == 1.2) 
      CI[1, 1] <- uniroot(tauRootDeltaInf, c(1e-4, max[1]), ss, ds,
                          pll.max, lim)$root
    else if(object$case == 1.12) {
      co[1] <- max[1]
      if(tauRoot(tau = 1e-4, ss, ds, sd, dd, pll.max, lim) < 0)
        ## Search for root below MLE
        CI[1, 1] <- uniroot(tauRoot, c(1e-4, co[1]), ss, ds, sd, dd,
                            pll.max, lim)$root
    }
    else if(is.na(co[2]) && object$case != 1.3) {
      if(co[1] == Inf)
        CI[1, 1] <- uniroot(tauRootDeltaInf, c(1e-4, max[1]), ss, ds,
                            pll.max, lim)$root
      else
        CI[1, 2] <- uniroot(tauRootDeltaInf, c(1e-4, max[1]), ss, ds,
                            pll.max, lim)$root
    }
    else {
      if(co[1] != Inf)
        ## Search for upper limit
        CI[1, 2] <- uniroot(tauRoot, c(co[1], max[1]), ss, ds, sd, dd,
                            pll.max, lim)$root
      else co[1] <- max[1]
      if(tauRoot(tau = 1e-4, ss, ds, sd, dd, pll.max, lim) < 0)
        ## Search for root below MLE
        CI[1, 1] <- uniroot(tauRoot, c(1e-4, co[1]), ss, ds, sd, dd,
                            pll.max, lim)$root
    }
  }
  if("delta" %in% parm){
    ## DELTA:
    if(co[2] != Inf)
      ## Search for upper limit
      CI[2, 2] <- uniroot(deltaRoot, c(co[2], max[2]), ss, ds, sd, dd,
                          pll.max, lim)$root
    else co[2] <- max[2]
    if(deltaRoot(delta = 1e-4, ss, ds, sd, dd, pll.max, lim) < 0)
      ## Search for root below MLE
      CI[2, 1] <- uniroot(deltaRoot, c(1e-4, co[2]), ss, ds, sd, dd,
                          pll.max, lim)$root
  }
  ## Return confidence interval:
  CI
}

contour.samediff <-
  function(x, norm = TRUE, level = c(.90, .95, .99), 
           lim = list(tau = c(1e-4, 10), delta = c(1e-4, 10),
             length = c(100, 100)), ...)
{
  da <- x$data
  dots <- list(...)
  Tau <- seq(lim$tau[1], lim$tau[2], length = lim$length[1])
  Delta <- seq(lim$delta[1], lim$delta[2], length = lim$length[2])
  ll <- outer(Tau, Delta, llSameDiff, ss=da[1], ds=da[2], sd=da[3],
              dd=da[4]) ## log-likelihood
  if(norm) { ## Normalized likelihood
    ll.max <- max(c(x$logLik), max(ll[!is.na(ll)])) 
    ll <- exp(ll - ll.max)
    contour(Tau, Delta, ll,
            xlab = expression(tau),
            ylab = expression(delta),
            levels = 1 - level,
            labels = level, ...)
  }
  else {
    contour(Tau, Delta, ll,
          xlab = expression(tau),
          ylab = expression(delta), ...)
  }
  invisible(list(tau = Tau, delta = Delta, ll = ll))
}
  
ell.delta <- function(tau.hat, delta, ss, ds, sd, dd) {
### Estimated log-likelihood for delta
  llSameDiff(tau.hat, delta, ss, ds, sd, dd)
}

ell.tau <- function(tau.hat, delta, ss, ds, sd, dd){
### Estimated log-likelihood for tau
  llSameDiff(tau.hat, delta, ss, ds, sd, dd)
}

sddvn <- 
  function(sn, n1, dn, n2)
{
#### d' and variance of d' for same-different method.
#### 
  s <- sn/n1
  d <- dn/n2
  p1 <- s + (1 - s)/2
  ta <- qnorm(p1) * sqrt(2)
  del <- function(delta, ta, d) 
    pnorm((ta - delta)/sqrt(2)) -
      pnorm(( - ta - delta)/sqrt(2)) - d
  delta <- uniroot(del, interval = c(0, 10), ta = ta, d = d) 
  delta <- delta[[1]]
  fi <- dnorm(qnorm(p1))
  fi <- fi^2
  btau <- (s * (1 - s))/(2 * n1 * fi)
  bvp2 <- (d * (1 - d))/n2
  a <- (ta - delta)/sqrt(2)
  b <- ( - ta - delta)/sqrt(2)
  dx <- - dnorm(a)/sqrt(2) + dnorm(b)/sqrt(2)
  dy <- dnorm(a)/sqrt(2) + dnorm(b)/sqrt(2)
  bdelta <- bvp2/dx^2 + ((dy^2) * btau)/(dx^2)
  dv <- c(delta, bdelta)
##  dv <- round(dv, 3)
  dv
}

plot.samediff <-
  function(x, main = TRUE, length = 1000,
           limits, fig = TRUE, ...)
{
  co <- coef(x)
  is.delta <- !is.na(co[2]) && co[2] != 0 && co[2] != Inf
  if(missing(limits)) {
    limits <- c(-4, 4)
    if(is.delta)
      limits[2] <- co[2] + 4
  }
  z <- seq(limits[1], limits[2], length.out = length)
  base.dist <- dnorm(z)
  object <- data.frame(z, base.dist)
  if(fig == TRUE) {
    main.txt <- ifelse(main,
                       paste("Distribution of sensory intensity for the",
                             "same-different test"), c("") )
    plot(z, base.dist, type="l", xlab = "Sensory Magnitude",
###       ylab = "Sensory Intensity",
         ylab = "", main = main.txt, las = 1, lty = 2, ...)
    if(is.delta) {
      object$delta.dist <- dnorm(z, mean = coef(x)[2])
      lines(z, object$delta.dist, col = "red", lty = 1, ...)
    }
  }
  invisible(object)
}

samediffSim <- function(n, tau, delta, Ns, Nd) {
  m <- match.call()
  m[[1]] <- as.name("list")
  eval.parent(m)
  pss <- 2 * pnorm(tau/sqrt(2)) - 1
  ss <- rbinom(n = n, size = Ns, prob = pss)
  ds <-  Ns - ss
  psd <- pnorm((tau - delta)/sqrt(2)) -
    pnorm((-tau - delta)/sqrt(2))
  sd <- rbinom(n = n, size = Nd, prob = psd)
  dd <- Nd - sd
  cbind(ss, ds, sd, dd)
}


samediffPwr <- function(n = 1000, tau, delta, Ns, Nd, alpha = 0.05) {
  m <- match.call()
  m[[1]] <- as.name("list")
  eval.parent(m)
  sds <- samediffSim(n, tau, delta, Ns, Nd)
  pval <- double(n)
  fun <- function(x) {
    sadi <- samediff(x[1], x[2], x[3], x[4])
    co <- coef(sadi)
    mll <- c(logLik(sadi))
    if(is.na(co[2]))
      pval.delta <- NA
    else if(co[2] == 0)
      pval.delta <- 1
    else  {
      ## If delta > 0 and potentially Inf:
      PllDelta <-
        pll.delta(, 1e-4, ss=x[1], ds=x[2], sd=x[3], dd=x[4])
      pval.delta <- 1 - pnorm(sqrt(2 * (mll - PllDelta)))
    }
    pval.delta
  }
  pval <- apply(sds, 1, fun)
  power <- mean(pval < alpha, na.rm = TRUE)
  power
}

