## standard double maximum functional

maxBB <- efpFunctional(lim.process = "Brownian bridge",
  functional = list(comp = function(x) max(abs(x)), time = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge", functional = "max", alt.boundary = FALSE, k = nproc))

maxBBI <- efpFunctional(lim.process = "Brownian bridge increments",
  functional = list(comp = function(x) max(abs(x)), time = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge increments", functional = "max", alt.boundary = FALSE, k = nproc))

maxBM <- efpFunctional(lim.process = "Brownian motion",
  functional = list(comp = function(x) max(abs(x)), time = max),
  boundary = function(x) 1 + 2*x,
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian motion", functional = "max", alt.boundary = FALSE, k = nproc))

maxBMI <- efpFunctional(lim.process = "Brownian motion increments",
  functional = list(comp = function(x) max(abs(x)), time = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian motion increments", functional = "max", alt.boundary = FALSE, k = nproc))


## range (over time) functional

rangeBB <- efpFunctional(lim.process = "Brownian bridge",
  functional = list(time = function(x) diff(range(x)), comp = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge", functional = "range", alt.boundary = FALSE, k = nproc))

rangeBBI <- efpFunctional(lim.process = "Brownian bridge increments",
  functional = list(time = function(x) diff(range(x)), comp = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge increments", functional = "range", alt.boundary = FALSE, k = nproc))

rangeBM <- efpFunctional(lim.process = "Brownian motion",
  functional = list(time = function(x) diff(range(x)), comp = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian motion", functional = "range", alt.boundary = FALSE, k = nproc))

rangeBMI <- efpFunctional(lim.process = "Brownian motion increments",
  functional = list(time = function(x) diff(range(x)), comp = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian motion increments", functional = "range", alt.boundary = FALSE, k = nproc))


## maximum L2-norm functional

maxL2BB <- efpFunctional(lim.process = "Brownian bridge",
  functional = list(comp = function(x) sum(x^2), time = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge", functional = "maxL2", alt.boundary = FALSE, k = nproc))


## Cramer-von Mises functional as employed in Nyblom-Hansen test

meanL2BB <- efpFunctional(lim.process = "Brownian bridge",
  functional = list(comp = function(x) sum(x^2), time = mean),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge", functional = "meanL2", alt.boundary = FALSE, k = nproc))


## Andrews' supLM test

supLM <- function(from = 0.15, to = NULL) {
  if(is.null(to)) to <- 1 - from
  stopifnot(from < 1, to < 1, from < to)
  lambda <- ((1-from)*to)/(from*(1-to))    
    
  wmax <- function(x) {
    n <- length(x)
    n1 <- floor(from * n)
    n2 <- floor(to * n)
    tt <- seq(along = x)/n
    x <- x[n1:n2]
    tt <- tt[n1:n2]
    x <- x/(tt * (1-tt))
    return(max(x))
  }

  computePval <- function(x, nproc = 1)
    pvalue.Fstats(x, type = "supF", k = nproc, lambda = lambda)
  
  computeCritval <- function(alpha, nproc = 1)
    uniroot(function(y) {computePval(y, nproc = nproc) - alpha}, c(0, 1000))$root

  boundary0 <- function(x) rep(1, length(x))

  plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
    xlab = NULL, ylab = NULL, main = x$type.name, ylim = NULL,
    boundary = TRUE, ...)
  {
    n <- x$nobs
    n1 <- floor(from * n)
    n2 <- floor(to * n)
    tt <- (1:n)/n

    bound <- computeCritval(alpha = alpha, nproc = NCOL(x$process)) * boundary0(0:n/n)
    if(is.null(xlab)) {
      if(!is.null(x$order.name)) xlab <- x$order.name
        else xlab <- "Time"
    }
    if(!identical(boundary, FALSE)) {
      if(isTRUE(boundary)) {
        bargs <- list(col = 2)
      } else if(!is.list(boundary)) {
        bargs <- list(col = boundary)
      } else {
        bargs <- boundary
      }
      boundary <- TRUE
    }

    if(aggregate) {
      proc <- zoo(rowSums(as.matrix(x$process)^2), time(x))
      proc <- proc[-1]
      bound <- zoo(bound, time(x)[-1])
      tt <- zoo(tt, time(x)[-1])
      
      proc <- proc[n1:n2,]
      bound <- bound[n1:n2,]
      tt <- as.numeric(tt[n1:n2,])
      proc <- proc/(tt * (1-tt))
      
      if(is.null(ylab)) ylab <- "Empirical fluctuation process"
      if(is.null(ylim)) ylim <- range(c(range(proc), range(bound)))
    
      plot(proc, xlab = xlab, ylab = ylab, main = main, ylim = ylim, ...)
      abline(0, 0)
      if(boundary) do.call("lines", c(list(bound), bargs))
    } else {
      if(is.null(ylim) & NCOL(x$process) < 2) ylim <- range(c(range(x$process), range(bound), range(-bound)))
      if(is.null(ylab) & NCOL(x$process) < 2) ylab <- "Empirical fluctuation process"

      panel <- function(x, ...)
      {
        lines(x, ...)
        abline(0, 0)
      }
      plot(x$process, xlab = xlab, ylab = ylab, main = main, panel = panel, ylim = ylim, ...)
    }
  }
  
  efpFunctional(lim.process = "Brownian bridge",
    functional = list(comp = function(x) sum(x^2), time = wmax),
    boundary = boundary0,
    computePval = computePval,
    computeCritval = computeCritval,
    plotProcess = plotProcess)
}

## (maximum) MOSUM functional

maxMOSUM <- function(width = 0.15) {

  mmax <- function(x) {
    n <- length(x)
    nh <- if(width < 1) floor(n * width) else round(width)
    x <- c(0, x)
    x <- x[-(1:nh)] - x[1:(n-nh+1)]
    max(abs(x))
  }

  computePval <- function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge increments", functional = "max",
      alt.boundary = FALSE, h = width, k = nproc)

  computeCritval <- function(alpha, nproc = 1)
    uniroot(function(y) {computePval(y, nproc = nproc) - alpha}, c(0, 100))$root

  boundary <- function(x) rep(1, length(x))

  plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
    xlab = NULL, ylab = NULL, main = x$type.name, ylim = NULL,
    boundary = TRUE, ...)
  {
    n <- x$nobs
    nh <- if(width < 1) floor(n * width) else round(width)
    nh2 <- floor(0.5 + nh/2)
    proc <- zoo(coredata(x$process[-(1:nh),]) - coredata(x$process[1:(n-nh+1),]),
      time(x$process)[nh2:(n - nh2 + 1)])

    if(!identical(boundary, FALSE)) {
      if(isTRUE(boundary)) {
        bargs <- list(col = 2)
      } else if(!is.list(boundary)) {
        bargs <- list(col = boundary)
      } else {
        bargs <- boundary
      }
      boundary <- TRUE
    }
    fboundary <- function(x) rep(1, length(x))
    bound <- computeCritval(alpha = alpha, nproc = NCOL(x$process)) * fboundary(nh2:(n - nh2 + 1)/n)
    bound <- zoo(bound, time(proc))

    if(is.null(xlab)) {
      if(!is.null(x$order.name)) xlab <- x$order.name
        else xlab <- "Time"
    }

    if(aggregate) {
      proc <- zoo(apply(as.matrix(proc), 1, function(x) max(abs(x))), time(proc))

      if(is.null(ylab)) ylab <- "Empirical fluctuation process"
      if(is.null(ylim)) ylim <- range(c(range(proc), range(bound)))
    
      plot(proc, xlab = xlab, ylab = ylab, main = main, ylim = ylim, ...)
      abline(0, 0)
      if(boundary) do.call("lines", c(list(bound), bargs))
    } else {
      if(is.null(ylim) & NCOL(x$process) < 2) ylim <- range(c(range(x$process), range(bound), range(-bound)))
      if(is.null(ylab) & NCOL(x$process) < 2) ylab <- "Empirical fluctuation process"

      panel <- function(x, ...)
      {
        lines(x, ...)
        abline(0, 0)
	if(boundary) {
	  do.call("lines", c(list(bound), bargs))
	  do.call("lines", c(list(-bound), bargs))
	}
      }
      plot(proc, xlab = xlab, ylab = ylab, main = main, panel = panel, ylim = ylim, ...)
    }
  }
  
  efpFunctional(lim.process = "Brownian bridge",
    functional = list(time = mmax, comp = max),
    boundary = boundary,
    computePval = computePval,
    computeCritval = computeCritval,
    plotProcess = plotProcess)
}

## efpFunctional() generator for categorical data
## (chi-squared statistic also used in MOB algorithm)

catL2BB <- function(freq)
{
  if(inherits(freq, "gefp")) freq <- as.factor(time(freq)[-1L])
  if(is.factor(freq)) freq <- prop.table(table(freq))
  freq <- freq / sum(freq)
  ncat <- length(freq)
  
  catdiffL2 <- function(x) {
    d <- diff(c(0, x[round(cumsum(freq) * length(x))]))
    sum(d^2 / freq)
  }
  
  computePval <- function(x, nproc = 1) pchisq(x, df = (ncat - 1) * nproc, lower.tail = FALSE)
  computeCritval <- function(alpha, nproc = 1) qchisq(alpha, df = (ncat - 1) * nproc, lower.tail = FALSE)

  efpFunctional(lim.process = "Brownian bridge",
    functional = list(time = catdiffL2, comp = sum),
    computePval = computePval,
    computeCritval = computeCritval)
}

## efpFunctional() generators for ordinal data
## ordL2BB: ordinal supLM-type statistic,
## ordwmax: weighted double-maximum type statistic

ordL2BB <- function(freq, nproc = NULL, nrep = 1e5, probs = c(0:84/100, 850:1000/1000), ...)
{
  if(inherits(freq, "gefp")) {
    if(is.null(nproc)) nproc <- NCOL(freq$process)
    freq <- prop.table(table(time(freq)[-1L]))
  } else if(is.factor(freq)) {
    if(is.null(nproc)) nproc <- 1L:20L
    freq <- prop.table(table(freq))
  } else {
    if(is.null(nproc)) nproc <- 1L:20L
  }
  freq <- freq / sum(freq)
  ncat <- length(freq)
  tcat <- cumsum(freq[-ncat])
  lab <- names(freq)
  if(is.null(lab)) lab <- as.character(1L:ncat)

  ## simulate quantiles of observations from asymptotic distribution
  probs <- probs[probs != 0]
  simordL2BB <- function(nproc, nrep, ...)
  {
    ## dimensions
    nbin <- ncat - 1L
    maxproc <- max(nproc)

    ## find multivariate normal covariance matrix between bins
    cmat <- matrix(NA, nbin, nbin)
    for(i in 1L:nbin) cmat[i:nbin, i] <- cmat[i, i:nbin] <- sqrt(tcat[i] * (1 - tcat[i:nbin]))/sqrt((1 - tcat[i]) * tcat[i:nbin])

    ## create block-diagonal matrix so we can simulate all parameters at once
    fullcmat <- matrix(0, nbin * maxproc, nbin * maxproc)
    for(j in 1L:maxproc) fullcmat[(nbin * (j - 1L) + 1L):(nbin * j), (nbin * (j - 1L) + 1L):(nbin * j)] <- cmat

    ## get samples from mvtnorm
    draws <- mvtnorm::rmvnorm(nrep, sigma = fullcmat, ...)

    ## now get L2 norm sequentially for each bin:
    ## columns 1,(nbin+1),(2nbin+1),...
    ## columns 2,(nbin+2),(2nbin+2),...
    colnums <- 1L:maxproc
    rval <- vector(mode = "list", length = nbin)
    for(j in 1L:nbin){
      cols <- (colnums - 1L) * nbin + j
      res <- t(apply(as.matrix(draws[, cols])^2, 1L, cumsum))
      if(maxproc == 1L) res <- matrix(res, nrep, 1L)
      rval[[j]] <- res[, nproc, drop = FALSE]
    }
    rval <- do.call("pmax", rval)
    rval <- rbind(0, apply(rval, 2, quantile, probs = probs))
    colnames(rval) <- nproc
    return(rval)
  }
  sim <- simordL2BB(nproc = nproc, nrep = nrep)

  ## compute p values and critical values from simulated quantiles
  computePval <- function(x, nproc = 1) {
    if(as.character(nproc) %in% colnames(sim)) {
      pfun <- approxfun(sim[, as.character(nproc)], 1 - c(0, probs))
      ifelse(x > max(sim[, as.character(nproc)]), 0, pfun(x))
    }
    else stop("insufficient simulated values: cannot compute p value")
  }
  computeCritval <- function(alpha, nproc = 1) {
    if(as.character(nproc) %in% colnames(sim)) {
      cfun <- approxfun(c(0, probs), sim[, as.character(nproc)])
      cfun(1 - alpha)
    } else stop("insufficient simulated values: cannot compute critical value")
  }

  ## "time" component of aggregation functional
  catwmax <- function(x) {
    n <- length(x)
    tt <- seq(along = x)/n
    ix <- round(tcat * n)
    x <- x[ix]
    tt <- tt[ix]
    x <- x/(tt * (1-tt))
    return(max(x))
  }

  ## constant boundaries
  boundary <- function(x) rep(1, length(x))

  ## visualization
  plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
    xlab = NULL, ylab = NULL, main = x$type.name, ylim = NULL,
    type = "b", boundary = TRUE, ...)
  {
    n <- x$nobs
    ix <- round(tcat * n)
    tt <- ix/n

    critval <- computeCritval(alpha = alpha, nproc = NCOL(x$process))
    bound <- critval * boundary(tt)
    if(!identical(boundary, FALSE)) {
      if(isTRUE(boundary)) {
        bargs <- list(col = 2)
      } else if(!is.list(boundary)) {
        bargs <- list(col = boundary)
      } else {
        bargs <- boundary
      }
      boundary <- TRUE
    }

    if(is.null(xlab)) {
      if(!is.null(x$order.name)) xlab <- x$order.name
        else xlab <- "Time"
    }

    if(aggregate) {
      proc <- zoo(rowSums(as.matrix(x$process)^2)[-1L][ix], time(x)[-1L][ix])
      bound <- zoo(bound, time(proc))
      tt <- zoo(tt, time(proc))      
      proc <- proc/(tt * (1-tt))
      
      if(is.null(ylab)) ylab <- "Empirical fluctuation process"
      if(is.null(ylim)) ylim <- range(c(range(proc), range(bound)))
    
      suppressWarnings(plot(proc, xlab = xlab, ylab = ylab, main = main, ylim = ylim,
        type = "b", axes = FALSE, ...))
      abline(0, 0)
      if(boundary) do.call("lines", c(list(bound), bargs))

      axis(1, at = unique(time(proc))[1L:(ncat - 1L)], labels = lab[1L:(ncat - 1L)])
      axis(2)
      box()
    } else {
      if(is.null(ylim) & NCOL(x$process) < 2L) ylim <- range(c(range(x$process)))
      if(is.null(ylab) & NCOL(x$process) < 2L) ylab <- "Empirical fluctuation process"

      panel <- function(x, ...)
      {
        lines(x, ...)
        abline(0, 0)
      }
      plot(x$process, xlab = xlab, ylab = ylab, main = main, panel = panel,
        ylim = ylim, ...)
    }
  }

  efpFunctional(lim.process = "Brownian bridge",
    functional = list(comp = function(x) sum(x^2), time = catwmax),
    boundary = boundary, plotProcess = plotProcess,
    computeCritval = computeCritval, computePval = computePval, ...)
}

ordwmax <- function(freq, algorithm = mvtnorm::GenzBretz(), ...)
{
  if(inherits(freq, "gefp")) {
    freq <- prop.table(table(time(freq)[-1L]))
  } else if(is.factor(freq)) {
    freq <- prop.table(table(freq))
  }
  freq <- freq / sum(freq)
  ncat <- length(freq)
  tcat <- cumsum(freq[-ncat])
  lab <- names(freq)
  if(is.null(lab)) lab <- as.character(1:ncat)

  catwmax <- function(x) {
    n <- length(x)
    tt <- seq(along = x)/n
    ix <- round(tcat * n)
    x <- x[ix]
    tt <- tt[ix]
    x <- x/sqrt(tt * (1-tt))
    return(max(abs(x)))
  }

  make_sigma <- function(tt) outer(tt, tt,
    function(x, y) sqrt(pmin(x, y) * (1 - pmax(x, y)) / ((pmax(x, y) * (1 - pmin(x, y))))))
  sigma <- make_sigma(tcat)
  
  computePval <- function(x, nproc = 1)
    (1 - mvtnorm::pmvnorm(lower = -x, upper = x, mean = rep(0, ncat - 1), sigma = sigma)^nproc)

  computeCritval <- function(alpha, nproc = 1)
    mvtnorm::qmvnorm((1 - alpha)^(1/nproc), tail = "both.tails", sigma = sigma)$quantile

  boundary <- function(x) rep(1, length(x))

  plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
    xlab = NULL, ylab = NULL, main = x$type.name, ylim = NULL,
    boundary = TRUE, type = "b", ...)
  {
    n <- x$nobs
    ix <- round(tcat * n)
    tt <- ix/n

    if(!identical(boundary, FALSE)) {
      if(isTRUE(boundary)) {
        bargs <- list(col = 2)
      } else if(!is.list(boundary)) {
        bargs <- list(col = boundary)
      } else {
        bargs <- boundary
      }
      boundary <- TRUE
    }
    cval <- if(boundary) computeCritval(alpha = alpha, nproc = NCOL(x$process)) else 0
    bound <- cval * boundary(tt)

    if(is.null(xlab)) {
      if(!is.null(x$order.name)) xlab <- x$order.name
        else xlab <- "Time"
    }

    proc <- zoo(as.matrix(x$process)[-1, , drop = FALSE][ix, , drop = FALSE], time(x)[-1][ix])
    bound <- zoo(bound, time(proc))
    tt <- zoo(tt, time(proc))	   
    proc <- proc/sqrt(tt * (1-tt))
    
    if(aggregate) proc <- zoo(apply(abs(proc), 1, max), time(proc))
    
    if(is.null(ylab)) ylab <- if(aggregate) "Empirical fluctuation process" else colnames(x$process)
    if(is.null(ylim)) ylim <- range(c(range(proc), 0, cval, if(aggregate) NULL else -cval))
    
    mypanel <- function(x, ...)
    {
      lines(x, ...)
      abline(0, 0)
      if(boundary) do.call("lines", c(list(bound), bargs))
      if(!aggregate & boundary) do.call("lines", c(list(-bound), bargs))
    }

    suppressWarnings(plot(proc, xlab = xlab, ylab = ylab, main = main, ylim = ylim, type = "b",
      panel = mypanel, axes = !aggregate, ...))

    if(aggregate) {
      axis(1, at = unique(time(proc))[1:(ncat - 1)], labels = lab[1:(ncat - 1)])
      axis(2)
      box()
    }
  }

  efpFunctional(lim.process = "Brownian bridge",
    functional = list(time = catwmax, comp = max),
    boundary = boundary,
    computePval = computePval,
    computeCritval = computeCritval,
    plotProcess = plotProcess, ...)
}
