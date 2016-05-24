### Check **ALL** implementations of the likelihood root statistic in
### profile functions etc. for correctness!

### In print.twoAC: Print the distribution of answers?


twoAC <-
  function(data, d.prime0 = 0, conf.level = 0.95,
           statistic = c("likelihood", "Wald"),
           alternative = c("two.sided", "less", "greater"),
           ## test = c("difference", "similarity"),
           ...)
### Have 'alternative' or 'test' argument?
### Include an "equiv[alent]" option to the alternative argument to
### allow for a TOST equivalence test?
### Gain argument method = c("preference", "discrimination") to
### disallow negative d-primes for discrimination (estimates, CIs,
### tests)?

### data - numeric vector of length 3 with non-negative entries.
###        Courced to integer.
### d.prime0 - numeric scalar
### conf.level - numeric scalar between 0 and 1.
{
  ## Initial argument matching:
  statName <- match.arg(statistic)
  alt <- match.arg(alternative)

  ## Initial argument consistency testing:
  if(!isTRUE(all.equal(round(data), data)))
    stop("non-integer numbers not allowed in 'data'")
  data <- as.integer(data)
  if(length(data) != 3)
    stop("'data' argument should be of length 3")
  stopifnot(is.numeric(conf.level) && length(conf.level) == 1 &&
            conf.level > 0 && conf.level < 1)
  stopifnot(is.numeric(d.prime0) && length(d.prime0) == 1)

  ## Get ML estimates, vcov, logLik:
  res <- estimate.2AC(data = data, vcov = TRUE, warn = FALSE)
  d.prime <- res$coefficients[2,1]
  se.d.prime <- res$coefficients[2,2]

  ## Add output to res:
  res <- c(list(alternative = alt, statistic = statName,
                conf.level = conf.level, call = match.call(),
                d.prime0 = d.prime0), res)
  class(res) <- "twoAC"

  ## Compute test statistic:
  if(statName == "likelihood")
    res$stat.value <-
      as.vector(LRtest.2AC(data, d.prime0 = d.prime0,
                           alternative = alt)[,"lroot"])
  if(statName == "Wald" && !is.null(res$vcov))
    res$stat.value <- (d.prime - d.prime0) / se.d.prime

  ## Compute p-value:
  if(!is.null(res$stat.value))
    res$p.value <- normalPvalue(statistic = res$stat.value,
                                alternative = alt)

  ## Get confidence intervals:
  if(!is.null(res$vcov))
    res$confint <- confint(res, parm = "d.prime", level = conf.level,
                           type = statName)

  ## return object:
  return(res)
}

print.twoAC <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
{
  ## Match test statistic:
  test.stat <- switch(x$statistic,
                      "likelihood" = "likelihood root statistic:",
                      "Wald" = "Wald statistic:")

  ## Print estimates:
  cat(paste("Results for the 2-AC protocol with data ",
            deparse(x$call$data), ":\n", sep = ""))
  ## printCoefmat(x$coefficients, tst.ind=integer(0), cs.ind=1:2,
  ##              digits=digits, P.values=FALSE, has.Pvalue=FALSE, ...)
### NOTE: Not using printCoefmat here because it wont print Inf and
### -Inf parameter estimates.
  print(x$coefficients, quote = FALSE, digits = digits, ...)

  ## Print confidence interval for d-prime:
  if(is.null(x$confint))
    cat("\nUse profile and confint methods to get confidence interval\n")
  else {
    cat(paste("\nTwo-sided ", round(100 * x$conf.level, 3),
              "% confidence interval for d-prime based on the\n", test.stat,
              "\n", sep = ""))
    ci <- x$confint
    colnames(ci) <- c("Lower", "Upper")
    print(ci, digits = digits, ...)
  }

  ## Print result of signifcance test:
  if(is.null(x$stat.value) && x$statistic == "Wald")
    cat("\nSignificance test not available - try with the likelihood statistic\n")
  ## this should never happen, but better safe than sorry:
  else if(is.null(x$stat.value) && x$statistic == "likelihood")
    cat("\nSignificance test not available\n")
  else {
    cat(paste("\nSignificance test:\n"))
    if(x$statistic == "Wald")
      cat(paste(" Wald statistic =", format(x$stat.value, digits),
                "p-value =", format.pval(x$p.value), "\n"))
    if(x$statistic == "likelihood")
      cat(paste(" Likelihood root statistic =",
                format(x$stat.value, digits),
                "p-value =", format.pval(x$p.value), "\n"))
    cat(" Alternative hypothesis: ")
    cat(paste("d-prime is", switch(x$alternative,
                                   "two.sided" = "different from",
                                   "less" = "less than",
                                   "greater" = "greater than"),
              format(x$d.prime0, digits), "\n"))
  }
  return(invisible(x))
}

profile.twoAC <-
  function(fitted, alpha = 1e-3, nSteps = 1e2, range, ...)
{
  ## Save d.prime for convenience:
  d.prime <- coef(fitted)[2, 1]
  if(is.na(d.prime))
    stop("profile is not available for d.prime = NA")
  ## Get range from Wald CI or over-rule with range argument:
  if(missing(range)) {
    ## Get range from Wald interval and supplied alpha:
    range <- as.vector(confint(fitted, parm = "d.prime",
                               level = 1 - alpha, type = "Wald"))
    if(!is.numeric(range) || !is.finite(range))
      stop(paste("could not determine 'range' from fitted object:",
                 "Please specify 'range'"))
  } else {
    stopifnot(is.numeric(range) && is.finite(range) &&
              length(range) >= 2)
    ## Warn if d.prime is finite and not in the range interval:
    if(is.finite(d.prime)) {
      if(d.prime >= max(range) || d.prime <= min(range))
        warning("d.prime should be in the interval given by range")
    }
  }
  dseq <- seq(from = min(range), to = max(range),
              length.out = nSteps)
  nll <- ## negative profile log-likelihood
    sapply(dseq, function(dd)
           optimize(nll.2AC, c(0, 10), d.prime = dd,
                    data = fitted$data)$objective)
  ## Compute the signed likelihood root statistic:
  Lroot <- sign(d.prime - dseq) * sqrt(2 * (fitted$logLik + nll))
  if(any(!is.finite(Lroot)))
    warning("invalid values of the likelihood profile occured")
  res <- data.frame("Lroot" = Lroot, "d.prime" = dseq)
  res <- res[order(res[,1]),]
  if(!all(diff(res[,2]) < 0))
    warning("likelihood is not monotonically decreasing from maximum,\n",
            "  so profile may be unreliable")
  prof <- vector("list", length = 1)
  names(prof) <- "d.prime"
  prof[[1]] <- res
  val <- structure(prof, original.fit = fitted)
  class(val) <- "profile.twoAC"
  return(val)
}

confint.twoAC <-
  function(object, parm, level = 0.95, type = c("likelihood", "Wald"),
           ...)
{
  type <- match.arg(type)
  if(type == "Wald") {
    cf <- coef(object)[,1]
    pnames <- names(cf)
    if (missing(parm))
        parm <- pnames
    else if (is.numeric(parm))
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(round(100*a, 1), "%")
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2L),
                dimnames = list(parm, pct))
    ses <- coef(object)[,2][parm]
    ci[] <- cf[parm] + ses %o% fac
  }  else
  if(type == "likelihood") {
    object <- profile(object, alpha = (1 - level) / 4, ...)
    ci <- confint(object, level=level, ...)
  }
  return(ci)
}

confint.profile.twoAC <-
  function(object, parm = "d.prime", level = 0.95, ...)
{
  if(parm != "d.prime")
    stop("Profile likelihood confidence interval only available for 'd.prime'")
  a <- (1-level)/2
  a <- c(a, 1-a)
  pct <- paste(round(100*a, 1), "%")
  ci <- array(NA, dim = c(1, 2),
              dimnames = list(parm, pct))
  cutoff <- qnorm(a)
  pro <- object[[ "d.prime" ]]
  sp <- spline(x = pro[, 2], y = pro[, 1])
  ci[1, ] <- approx(-sp$y, sp$x, xout = cutoff)$y
  return(ci)
}

plot.profile.twoAC <-
  function(x, level = c(0.95, 0.99), Log = FALSE, relative = TRUE,
           fig = TRUE, n = 1e3, ..., ylim = NULL)
{
  ML <- attr(x, "original.fit")$logLik
  lim <- sapply(level, function(x)
                exp(-qchisq(x, df = 1)/2) )
  pro <- x[[ "d.prime" ]]
  sp <- spline(x = pro[, 2], y = pro[, 1], n = n)
  sp$y <- -sp$y^2/2
  if(relative && !Log) {
    sp$y <- exp(sp$y)
    ylab <- "Relative profile likelihood"
    dots <- list(...)
    if(missing(ylim))
      ylim <- c(0, 1)
  }
  if(relative && Log) {
    ylab <- "Relative profile log-likelihood"
    lim <- log(lim)
  }
  if(!relative & Log) {
    sp$y <- sp$y + ML
    ylab <- "profile Log-likelihood"
    lim <- ML + log(lim)
  }
  if(!relative & !Log) {
    stop("Not supported: at least one of 'Log' and 'relative' ",
         "have to be TRUE")
    sp$y <- exp(sp$y + ML)
    ylab <- "Profile likelihood"
    lim <- exp(ML + log(lim))
  }
  x[[ "d.prime" ]] <- sp
  if(fig) {
    plot(sp$x, sp$y, type = "l", ylim = ylim,
         xlab = "d.prime", ylab = ylab, ...)
    abline(h = lim)
  }
  attr(x, "limits") <- lim
  return(invisible(x))
}

nll.2AC <- function(tau, d.prime = 0, data)
### Computes the negative log-likelihood of the 2-AC protocol
{
  ## get probability vector (prob) from 2-AC parameters:
  p1 <- pnorm(-tau, d.prime, sqrt(2))
  p12 <- pnorm(tau, d.prime, sqrt(2))
  prob <- c(p1, p12 - p1, 1 - p12)
  prob[prob <= 0] <- 1 ## to evaluate log safely

  ## Evaluate negative log likelihood:
  -sum(data * log(prob))
}

estimate.2AC <- function(data, vcov = TRUE, warn = TRUE)
### Estimate the parameters, vcov and standard errors of the 2-AC
### model given a 3-vector of data.

### data - vector of data of length 3
### vcov - should the variance-covariance matrix, and hence the
### standard errors of the parameters be computed?
{
  ## test data argument?
  vcov <- as.logical(vcov)
  stopifnot(is.logical(vcov))

  ## Define negative log-likelihood:
  nll <- function(par) {
    tau <- par[1]
    d.prime <- par[2]
    p1 <- pnorm(-tau, d.prime, sqrt(2))
    p12 <- pnorm(tau, d.prime, sqrt(2))
    prob <- c(p1, p12 - p1, 1 - p12)
    prob[prob <= 0] <- 1 ## to evaluate log safely
    ## Evaluate negative log likelihood:
    -sum(data * log(prob))
  }

  ## Get ML estimates:
  x <- data
  if(x[1] > 0 && x[2] == 0 && x[3] == 0) { # case 1
    tau <- 0
    d.prime <- -Inf
  }
  else if(x[1] == 0 && x[2] > 0 && x[3] == 0) { # case 2
    tau <- NA
    d.prime <- NA
  }
  else if(x[1] == 0 && x[2] == 0 && x[3] > 0) { # case 3
    tau <- 0
    d.prime <- Inf
  }
  else if(x[1] > 0 && x[2] > 0 && x[3] == 0) { # case 4
    d.prime <- -Inf
    tau <- NA
  }
  else if(x[1] == 0 && x[2] > 0 && x[3] > 0) { # case 5
    d.prime <- Inf
    tau <- NA
  }
  else { # case 0 and 6
    prob <- data / sum(data)
    gamma <- cumsum(prob)[-3]
    z <- qnorm(gamma) * sqrt(2)
    tau <- (z[2] - z[1]) / 2
    d.prime <- -z[1] - tau
  }

  ## Get the log likelihood at the MLE:
  prob <- data / sum(data)
  prob[prob <= 0] <- 1 ## to evaluate log safely
  logLikMax <- sum(data * log(prob))

  ## Save estimates in coefficient table
  coef <- matrix(NA, 2, 2, dimnames = list(c("tau", "d.prime"),
                            c("Estimate", "Std. Error")))
  coef[,1] <- c(tau, d.prime)
  res <- list(coefficients = coef, data = data,
              logLik = logLikMax)

  ## Get Hessian, vcov and standard errors:
  if(vcov) {
    makeWarn <- TRUE
    ## If all coef are finite and tau < 0:
    if(all(is.finite(coef[,1])) && tau > 0) {
      makeWarn <- FALSE
      hess <- hessian(nll, x = c(tau, d.prime), method = "Richardson",
                      method.args = list())
      vcov <- try(solve(hess), silent = TRUE)
      ## If hess is not invertible:
      if(class(vcov) != "try-error") {
        makeWarn <- TRUE
        res$coefficients[,2] <- sqrt(diag(vcov))
        res$vcov <- vcov
      }
    }
    ## Not all coef are finite or tau <= 0:
    if(warn && makeWarn)
      warning("vcov and standard errors are not available",
              call. = FALSE)
  } ## end vcov

  ## Return results
  return(res)
}

vcov.twoAC <- function(object, ...) object$vcov
logLik.twoAC <- function(object, ...) object$logLik

twoACpwr <-
  function(tau, d.prime, size, d.prime0 = 0, alpha = 0.05, tol = 1e-5,
           return.dist = FALSE, statistic = "likelihood",
           alternative = c("two.sided", "less", "greater"))
  ## allow for statistic = c("likelihood", "Wald")?
{
  ## Match arguments:
  alternative <- match.arg(alternative)
  statistic <- match.arg(statistic)
  return.dist <- as.logical(return.dist)

  ## Initial testing:
  ## test for type/mode, length, valid values, integer values
  stopifnot(is.numeric(d.prime) && length(d.prime) == 1)
  stopifnot(is.numeric(tau) && length(tau) == 1 && tau > 0)
  stopifnot(is.numeric(d.prime0) && length(d.prime0) == 1)
  stopifnot(is.numeric(size) &&
            isTRUE(all.equal(round(size), size)) &&
            size < 5e3 && size >= 2)
  stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
            alpha > 0 && alpha < 1)
  stopifnot(is.numeric(tol) && length(tol) == 1,
            tol >= 0 && tol < 1)
  stopifnot(is.logical(return.dist))

  ## All possible samples:
  n <- as.integer(size)
  Y <- do.call(rbind, lapply(0:n, function(j) cbind(n-j, 0:j, j:0)))

  ## Get pi-vector from tau and d.prime:
  p1 <- pnorm(-tau, d.prime, sqrt(2))
  p2 <- pnorm(tau, d.prime, sqrt(2)) - p1
  p3 <- 1 - p1 - p2
  pvec <- c(p1, p2, p3)

  ## Distribution of samples and p-values:
  dmult <- apply(Y, 1, function(x) dmultinom(x, prob = pvec))

  ## Discard tail of distribution - keep only those obs with large
  ## enough probability mass to avoid computing p-values for samples
  ## that never occur anyway:
  sdmult <- sort(dmult)
  no.discard <- sum(cumsum(sdmult) < tol) ## no. obs to discard
  keep <- dmult > sdmult[no.discard] # obs. to keep indicator
  if(no.discard == 0) keep <- rep(TRUE, length(dmult))
  dmult <- dmult[keep]

  ## Compute p-values:
  pvals <- LRtest.2AC(Y[keep, ], d.prime0 = d.prime0,
                      alternative = alternative)[,1]
  ## pvals <- apply(Y[keep, ], 1, function(x) LRtest.2AC(x)[2])

  ## Compute the cumulative distribution of p-values:
  df <- data.frame(dens=dmult, p.value=pvals, Y=Y[keep, ])
  df <- df[order(pvals),]
  dist <- cumsum(df[,1])
  pvals <- df[,2]
  if(return.dist)
    return(structure(data.frame(dist=dist, df), row.names=1:nrow(df)))
  if(any(pvals <= alpha)) {
    power <- max(dist[pvals <= alpha])
    actual.alpha <- max(pvals[pvals <= alpha])
  }
  else {
    power <- 0
    actual.alpha <- 0
  }
  return(data.frame("power" = power, "actual.alpha" = actual.alpha,
                    "samples" = nrow(Y), "discarded" = no.discard,
                    "kept" = nrow(Y) - no.discard,
                    "p" = round(matrix(pvec, nrow = 1), 4)))
}

LRtest.2AC <-
  function(data, d.prime0 = 0,
           alternative = c("two.sided", "less", "greater")
           ## statistic = c("likelihood", "Wald")
           )
### Computes the signed likelihood root statistic and p-value for the
### indicated significance test. Care is taken that the correct
### statistic is computed in boundary cases with null cells.

### data - a 3-vector or a 3-col matrix where each row is taken as a
### trial
### d.prime0 - value of d.prime under the null hypothesis, scalar.
{
  ## Argument matching:
  alternative <- match.arg(alternative)

  ## Initial testing:
  stopifnot(is.numeric(d.prime0) && length(d.prime0) == 1)
  stopifnot(is.numeric(data) && all(data >= 0))
  if(is.null(dim(data)) && length(data) == 3) dim(data) <- c(1, 3)
  if(is.data.frame(data)) data <- as.matrix(data)
  if(!(is.matrix(data) &&  ncol(data) == 3))
    stop("'data' should be a vector of length 3 or a matrix with 3 cols")
  ## test that isTRUE(all.equal(round(data), data)) ignored here
  ## because of speed.

  ## Make indicator variable with the seven possible data cases:
  cases <- integer(length = nrow(data))
  cases[data[,3] == 0] <- 4L
  cases[data[,1] == 0] <- 5L
  cases[data[,2] == 0] <- 6L
  cases[data[,2] == 0 & data[,3] == 0] <- 1L
  cases[data[,1] == 0 & data[,3] == 0] <- 2L
  cases[data[,1] == 0 & data[,2] == 0] <- 3L

  ## Compute the signed likelihood root statistic in all seven cases:
  lroot <- numeric(length = nrow(data))
  if(sum(ind <- (cases == 1)))
    lroot[ind] <- sapply(which(ind), function(i) {
      nll0 <- nll.2AC(0, d.prime = d.prime0, data = data[i,])
      -sqrt(2 * nll0) ## lroot statistic
    })
  ## Not needed since lroot = 0 by default:
  ## if(sum(ind <- (cases == 2)))
  ##   lroot[ind] <- 0
  if(sum(ind <- (cases == 3)))
    lroot[ind] <- sapply(which(ind), function(i) {
      nll0 <- nll.2AC(0, d.prime = d.prime0, data = data[i,])
      sqrt(2 * nll0) ## lroot statistic
    })
  if(sum(ind <- (cases == 4)))
    lroot[ind] <- sapply(which(ind), function(i) {
      nll0 <- optimize(nll.2AC, c(0, 10), d.prime = d.prime0,
                       data = data[i, ])$objective
      prob <- data[i, ] / sum(data[i, ])
      prob[prob <= 0] <- 1 ## to evaluate log safely
      llMax <- sum(data[i, ] * log(prob))
      -sqrt(2 * (llMax + nll0)) ## lroot statistic
    })
  if(sum(ind <- (cases == 5)))
    lroot[ind] <- sapply(which(ind), function(i) {
      nll0 <- optimize(nll.2AC, c(0, 10), d.prime = d.prime0,
                       data = data[i, ])$objective
      prob <- data[i, ] / sum(data[i, ])
      prob[prob <= 0] <- 1 ## to evaluate log safely
      llMax <- sum(data[i, ] * log(prob))
      sqrt(2 * (llMax + nll0)) ## lroot statistic
    })
  if(sum(ind <- (cases %in% c(0, 6))))
    lroot[ind] <- sapply(which(ind), function(i) {
      nll0 <- optimize(nll.2AC, c(0, 10), d.prime = d.prime0,
                       data = data[i, ])$objective
      prob <- data[i, ] / sum(data[i, ])
      d.prime <- -sum(qnorm(cumsum(prob)[-3]) / sqrt(2))
      prob[prob <= 0] <- 1 ## to evaluate log safely
      llMax <- sum(data[i, ] * log(prob))
      sign(d.prime - d.prime0) * sqrt(2 * (llMax + nll0)) ## lroot statistic
    })

  ## Compute p-value given the alternative hypothesis:
  p.value <-
    switch(alternative,
           "greater" = pnorm(lroot, lower.tail = FALSE),
           "less" = pnorm(lroot, lower.tail = TRUE),
           "two.sided" = 2 * pnorm(abs(lroot), lower.tail = FALSE))

  return(cbind(p.value = p.value, lroot = lroot))
}

pearsonPwr <-
  function(tau, d.prime, size, tol=1e-5, return.dist=FALSE, alpha=0.05,
           alternative = c("two.sided", "less", "greater"))

### In this function d.prime0, (i.e., d.prime under the null
### hypothesis) is fixed at zero
{
  alternative <- match.arg(alternative)
  ## All possible samples:
  n <- as.integer(size)
  Y <- do.call(rbind, lapply(0:n, function(j) cbind(n-j, 0:j, j:0)))

  ## Get pi-vector from tau and d.prime:
  p1 <- pnorm(-tau, d.prime, sqrt(2))
  p2 <- pnorm(tau, d.prime, sqrt(2)) - p1
  p3 <- 1 - p1 - p2
  pvec <- c(p1, p2, p3)

  ## Distribution of samples and p-values:
  dmult <- apply(Y, 1, function(x) dmultinom(x, prob = pvec))
  sdmult <- sort(dmult)
  no.discard <- sum(cumsum(sdmult) < tol) ## no. obs to discard
  keep <- dmult > sdmult[no.discard] # obs. to keep indicator
  if(no.discard == 0) keep <- rep(TRUE, length(dmult))
  dmult <- dmult[keep]

  ## pvals <- apply(Y[keep, ], 1, Pears)
  ## pvals[is.na(pvals)] <- 1
  pvals <- apply(Y[keep, ], 1, function(x) Potter(x, alternative=alternative))

  df <- data.frame(dens=dmult, p.value=pvals, Y=Y[keep, ])
  df <- df[order(pvals),]
  dist <- cumsum(df[,1])
  pvals <- df[,2]
  if(return.dist)
    return(structure(data.frame(dist=dist, df), row.names=1:nrow(df)))
  if(any(pvals <= alpha)) {
    power <- max(dist[pvals <= alpha])
    actual.alpha <- max(pvals[pvals <= alpha])
  }
  else {
    power <- 0
    actual.alpha <- 0
  }
  return(data.frame("power" = power, "actual.alpha" = actual.alpha,
                    "samples" = nrow(Y), "discarded" = no.discard,
                    "kept" = nrow(Y) - no.discard,
                    "p" = round(matrix(pvec, nrow = 1), 4)))
}

Pears <- function(x) {
  y <- (x[1] + x[3])/2
  X <- sum((x[-2] - y)^2 / y)
  pchisq(X, df=1, lower.tail=FALSE)
}

Potter <- function(x, alternative) {
  n <- sum(x[c(1,3)])
  if(n == 0) return(1)
  X <- (x[1] - x[3]) / sqrt(n)
  p.value <-
    switch(alternative,
           "greater" = pnorm(X, lower.tail = FALSE),
           "less" = pnorm(X, lower.tail = TRUE),
           "two.sided" = 2 * pnorm(abs(X), lower.tail = FALSE))
  p.value
}


binPwr <-
  function(tau, d.prime, size, tol=1e-5, return.dist=FALSE, alpha=0.05,
           alternative = c("two.sided", "less", "greater"))
### In this function d.prime0, (i.e., d.prime under the null
### hypothesis) is fixed at zero
{
  alternative <- match.arg(alternative)
  ## All possible samples:
  n <- as.integer(size)
  Y <- do.call(rbind, lapply(0:n, function(j) cbind(n-j, 0:j, j:0)))

  ## Get pi-vector from tau and d.prime:
  p1 <- pnorm(-tau, d.prime, sqrt(2))
  p2 <- pnorm(tau, d.prime, sqrt(2)) - p1
  p3 <- 1 - p1 - p2
  pvec <- c(p1, p2, p3)

  ## Distribution of samples and p-values:
  dmult <- apply(Y, 1, function(x) dmultinom(x, prob = pvec))
  sdmult <- sort(dmult)
  no.discard <- sum(cumsum(sdmult) < tol) ## no. obs to discard
  keep <- dmult > sdmult[no.discard] # obs. to keep indicator
  if(no.discard == 0) keep <- rep(TRUE, length(dmult))
  dmult <- dmult[keep]

  pvals <- apply(Y[keep, ], 1, function(x) binExact(x, alternative=alternative))
  pvals[is.na(pvals)] <- 1

  df <- data.frame(dens=dmult, p.value=pvals, Y=Y[keep, ])
  df <- df[order(pvals),]
  dist <- cumsum(df[,1])
  pvals <- df[,2]
  if(return.dist)
    return(structure(data.frame(dist=dist, df), row.names=1:nrow(df)))
  if(any(pvals <= alpha)) {
    power <- max(dist[pvals <= alpha])
    actual.alpha <- max(pvals[pvals <= alpha])
  }
  else {
    power <- 0
    actual.alpha <- 0
  }
  return(data.frame("power" = power, "actual.alpha" = actual.alpha,
                    "samples" = nrow(Y), "discarded" = no.discard,
                    "kept" = nrow(Y) - no.discard,
                    "p" = round(matrix(pvec, nrow = 1), 4)))
}

binExact <- function(x, alternative) {
  n <- x[1] + x[3]
  ifelse(n == 0, 1, binom.test(x[1], x[1]+x[3], alternative=alternative)$p.value)
}

clm2twoAC <- function(object, ...) {
  if(inherits(object, c("clm", "clmm")))
    theta <- object$alpha
  else if(inherits(object, c("clm2", "clmm2")))
    theta <- object$Theta
  else
    stop("'object' not of appropriate class")
  tab <- coef(summary(object))
  stopifnot(length(theta) == 2)
  mat <- array(NA, dim = c(2, 4))
  rownames(mat) <- c("tau", "d-prime")
  colnames(mat) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  ## parameter estimates:
  mat[1, 1] <- (theta[2] - theta[1]) / sqrt(2)
  mat[2, 1] <- (-theta[2] - theta[1]) / sqrt(2)
  VCOV <- vcov(object)[1:2, 1:2]
  ## standard errors:
  mat[1, 2] <- sqrt((VCOV[1, 1] + VCOV[2,2] - 2* VCOV[2, 1])/2)
  mat[2, 2] <- sqrt((VCOV[1, 1] + VCOV[2,2] + 2* VCOV[2, 1])/2)
  ## z-values and p-values:
  mat[,3] <- mat[,1] / mat[,2]
  mat[,4] <- 2*pnorm(abs(mat[, 3]), lower.tail=FALSE)
  ## add additional rows to coefficient matrix
  if(dim(tab)[1] > 2) {
    tmp <- tab[-(1:2), , drop = FALSE]
    ## scale estimates and standard errors with sqrt(2):
    tmp[,1:2] <- tmp[,1:2] * sqrt(2)
    mat <- rbind(mat, tmp)
    rownames(mat)[-(1:2)] <- rownames(tab)[-(1:2)]
  }
  mat <- as.data.frame(mat)
  mat[,4] <- format.pval(mat[,4])
  return(mat)
}

exact.2AC <-
  function(x, d.prime0=0, tau0=NULL, tol=1e-5,
           alternative = c("two.sided", "less", "greater"), ...)
{
  alternative <- match.arg(alternative)
  coefs <- coef(estimate.2AC(x, vcov=FALSE))[,1]
  dHat <- coefs[2]
  n <- as.integer(sum(x))
  Y <- do.call(rbind, lapply(0:n, function(j) cbind(n-j, 0:j, j:0)))

  ## get pvec under the null hypothesis:
  if(is.null(tau0)) ## tau0 <- coefs[1]
    tau0 <- optimize(nll.2AC, interval=c(0, 10),
                     d.prime=d.prime0, data=x)$minimum

  ## Get pi-vector from tau and d.prime:
  p1 <- pnorm(-tau0, d.prime0, sqrt(2))
  p2 <- pnorm(tau0, d.prime0, sqrt(2)) - p1
  p3 <- 1 - p1 - p2
  pvec <- c(p1, p2, p3)

  ## Distribution of samples under the null:
  dmult <- apply(Y, 1, function(x) dmultinom(x, prob = pvec))

  ## Discard tail of distribution - keep only those obs with large
  ## enough probability mass to avoid computing d.prime for samples
  ## that never occur anyway:
  sdmult <- sort(dmult)
  no.discard <- sum(cumsum(sdmult) < tol) ## no. obs to discard
  keep <- dmult > sdmult[no.discard] # obs. to keep indicator
  if(no.discard == 0) keep <- rep(TRUE, length(dmult))
  dmult <- dmult[keep]

  ## d.primes for each possible outcome:
  d.primes <- apply(Y[keep, ], 1, function(x)
                    coef(estimate.2AC(x, vcov=FALSE))[2,1] )
  df <- data.frame(dens=dmult, d.prime=d.primes)
  df <- df[!is.na(d.primes),]

  p.value <-
    switch(alternative,
           "greater" = with(df, sum(dens[d.prime >= dHat])),
           "less" = with(df, sum(dens[d.prime <= dHat])),
           "two.sided" = with(df,
             sum(dens[d.prime <= -dHat | d.prime >= dHat])))
### FIXME: should be abs(dHat) here?
  res <-
    data.frame("p.value"=p.value, "d.prime0"=d.prime0,
               "tau0"=tau0, "d.prime.hat"=dHat, "samples"=nrow(Y),
               "discarded" = no.discard, "kept"=nrow(Y) - no.discard)
  row.names(res) <- ""
  res
}



##  LRtest.2AC.old <- function(x, ...) {
##  ### With the fomulation below, this function is only applicable when
##  ### d.prime0 = 0, i.e. when value of d.prime under the null hypothesis
##  ### is zero.
##  ### x: vector of length 3 containing the data for the 2-AC protocol.
##  ### d.prime0: Value under the null hypothesis
##  ### alternative = c("two-sided", "greater", "less")
##  ### statistic = c("likelihood", "Wald")
##    nll.tau <- function(tau, d = 0, data) {
##      p1 <- pnorm(-tau, d, sqrt(2))
##      p12 <- pnorm(tau, d, sqrt(2))
##      prob <- c(p1, p12 - p1, 1 - p12)
##      prob[prob <= 0] <- 1 ## to evaluate log safely
##      -sum(data * log(prob)) # negative log likelihood
##    }
##
##    ## log-likelihood at ML estimates:
##    prob <- x / sum(x)
##    prob[prob <= 0] <- 1 ## to evaluate log safely
##    llMax <- sum(x * log(prob))
##
##    ## log-likelihood and tau at d-prime = 0:
##    if(x[2] == 0) { # & other x arbitrary
##      tau <- 0
##      nll.d0 <- nll.tau(0, 0, x)
##    }
##    else if(x[1] == 0 && x[3] == 0) { # & x[2] > 0
##      tau <- Inf
##      nll.d0 <- 0
##    }
##    else {
##      fit <- optimize(nll.tau, c(0, 10), data = x)
##      nll.d0 <- fit$objective
##      tau <- fit$minimum
##    }
##    ## Signed likelihood root statistic:
##    ## lroot <- sign(d.prime.hat - d.prime0) * sqrt(2*(llMax + nll.d0))
##    lroot <- abs(sqrt(2*(llMax + nll.d0))) ## absolute likelihood root
##    ## statistic
##    ## if(alternative == "two-sided")
##    p.value <- 2 * pnorm(abs(lroot), lower.tail = FALSE)
##    ## if(alternative == "difference")
##    ##   p.value <- pnorm(lroot, lower.tail = FALSE)
##    ## if(alternative == "similarity")
##    ##   p.value <- pnorm(lroot, lower.tail = TRUE)
##    return(c(lroot, p.value, tau))
##  }

## getAnywhere("coef.default")
## getAnywhere("logLik.default")
## getAnywhere("AIC.default")
## getAnywhere("BIC.default")
## getAnywhere("logLik")
