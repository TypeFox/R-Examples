gefp <- function(...,
  fit = glm, scores = estfun, vcov = NULL,
  decorrelate = TRUE, sandwich = TRUE,
  order.by = NULL, fitArgs = NULL, parm = NULL, data = list())
{
  ## avoid name clashes of vcov argument and vcov() function
  vcov. <- vcov

  if(is.null(fit)) {
    fm <- list(...)
    if(length(fm) > 1) warning("more than one argument supplied in `...', only the first is used")
    fm <- fm[[1]]
  } else {
    if(is.null(fitArgs)) fm <- fit(..., data = data)	
      else fm <- do.call("fit", c(..., fitArgs, list(data = data)))
  }
    
  psi <- scores(fm)
  n <- NROW(psi)
  k <- NCOL(psi)

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.frame(order.by, data = data)
      z <- as.numeric(as.vector(z[, ncol(z)]))
      order.name <- deparse(order.by[[2]])
    } else {
      z <- order.by
      order.name <- deparse(substitute(order.by))
    }
    index <- order(z)
    psi <- psi[index, , drop = FALSE]
    z <- z[index]
    if(is.factor(z)) z <- as.numeric(z) ## FIXME: deal better with factor orderings?
  } else {
    index <- 1:n
    if(is.ts(psi)) z <- time(psi)
      else if(is.zoo(psi)) z <- time(psi)
      else if(is.ts(data)) z <- time(data)
      else if(is.zoo(data)) z <- time(data)
      else ## support for time series regression in lm()
           ## has been limited since R 2.0.0. Until there is
           ## a better option, there is an ugly hack here,
           ## for `guessing' ts attributes.
	  if(any(c("lm", "glm", "rlm") %in% class(fm)) && !is.null(formula(fm))) {
            form <- formula(fm)
            env <- environment(form)
	    if(is.null(fit)) {
	      dat <- fm$call[["data"]]
	      if(!is.null(dat)) dat <- get(deparse(dat), envir = env)
	    } else dat <- NULL
	    if(is.ts(dat)) {
	      z <- time(dat)
	    } else {
  	      if(missing(data)) data <- if(is.null(dat)) env else dat
              orig.y <- eval(attr(terms(form), "variables")[[2]], data, env)
              if(is.ts(orig.y)) z <- time(orig.y)
	        else z <- index/n
	    }	    
          } else z <- index/n
    order.name <- "Time"
  }
  
  psi <- as.matrix(psi)

  if(inherits(z, "POSIXt"))
    z <- suppressWarnings(c(z[1] + as.numeric(difftime(z[1], z[2], units = "secs")), z))
  else
    z <- suppressWarnings(c(z[1] - as.numeric(diff(z[1:2])), z))

  process <- psi/sqrt(n)

  if(is.null(vcov.))
    J12 <- root.matrix(crossprod(process))
  else {
    if(sandwich) {
      Q <- chol2inv(chol(bread(fm)/n))
      J12 <- root.matrix((Q %*% vcov.(fm, order.by = order.by, data = data) %*% Q)/n)
    } else {
      J12 <- root.matrix(vcov.(fm, order.by = order.by, data = data))
    }
  }

  process <- rbind(0, process)
  process <- apply(process, 2, cumsum)

  if(decorrelate) process <- t(chol2inv(chol(J12)) %*% t(process))
    else {
      process <- t(1/diag(J12) * t(process))
      if(length(parm) > 1) stop("limiting process is not a Brownian bridge")
    }

  colnames(process) <- colnames(psi)
  if(!is.null(parm)) process <- process[, parm]

  retval <- list(process = suppressWarnings(zoo(process, z)),
                 nreg = k,
                 nobs = n,
                 call = match.call(),
		 fit = fit,
		 scores = scores,
		 fitted.model = fm,
                 par = NULL,
                 lim.process = "Brownian bridge",
		 type.name = "M-fluctuation test",
		 order.name = order.name,
                 J12 = J12)

  class(retval) <- "gefp"
  return(retval)
}

plot.gefp <- function(x, alpha = 0.05, functional = maxBB, ...)
{
  if(is.null(functional)) functional <- "max"
  if(is.character(functional)) {
    functional <- match.arg(functional, c("max", "range", "maxL2", "meanL2"))
    lim.process <- switch(x$lim.process,
      "Brownian motion" = "BM",
      "Brownian motion increments" = "BMI",
      "Brownian bridge" = "BB",
      "Brownian bridge increments" = "BBI")
    functional <- get(paste(functional, lim.process, sep = ""), pos = "package:strucchange")
  }
  if(!("efpFunctional" %in% class(functional)))
    stop(paste(dQuote("functional"), "has to be of class", sQuote("efpFunctional")))
  if(functional$lim.process != x$lim.process)
    stop(paste("limiting process of", dQuote("functional"), "does not match that of", dQuote("x")))
  functional$plotProcess(x, alpha = alpha, ...)
}

sctest.gefp <- function(x, functional = maxBB, ...)
{
  if(is.character(functional)) {
    functional <- match.arg(functional, c("max", "range", "maxL2", "meanL2"))
    lim.process <- switch(x$lim.process,
      "Brownian motion" = "BM",
      "Brownian motion increments" = "BMI",
      "Brownian bridge" = "BB",
      "Brownian bridge increments" = "BBI")
    functional <- get(paste(functional, lim.process, sep = ""), pos = "package:strucchange")
  }
  if(!("efpFunctional" %in% class(functional)))
    stop(paste(dQuote("functional"), "has to be of class", sQuote("efpFunctional")))
  if(functional$lim.process != x$lim.process)
    stop(paste("limiting process of", dQuote("functional"), "does not match that of", dQuote("x")))
  stat <- functional$computeStatistic(x$process)
  pval <- functional$computePval(stat, NCOL(x$process))
  names(stat) <- "f(efp)"
  rval <- list(statistic = stat,
               p.value = pval,
	       method = x$type.name,
	       data.name = deparse(substitute(x)))
  class(rval) <- "htest"
  return(rval)
}

sctest.default <- function(x, order.by = NULL, functional = maxBB,
  vcov = NULL, scores = estfun, decorrelate = TRUE, sandwich = TRUE, parm = NULL,
  plot = FALSE, from = 0.1, to = NULL, nobs = NULL, nrep = 50000, width = 0.15,
  xlab = NULL, ...)
{
  ## object names
  nam <- deparse(substitute(x))
  if(is.null(xlab)) {
    if(is.null(order.by)) {
      xlab <- "Time"
    } else {
      xlab <- deparse(substitute(order.by))
      xlab <- tail(unlist(strsplit(xlab, "$", fixed = TRUE)), 1)
    }
  }
  
  ## convenience option to employ information matrix
  ## (sensible for maximum likelihood fits only, otherwise needs scaling)
  vcov. <- vcov
  if(identical(vcov., "info")) {
    vcov0 <- if("stats4" %in% loadedNamespaces()) stats4::vcov else stats::vcov
    nobs0   <- function(x, ...) {
      nobs1 <- if("stats4" %in% loadedNamespaces()) stats4::nobs else stats::nobs
      nobs2 <- function(x, ...) NROW(residuals(x, ...))
      rval <- try(nobs1(x, ...), silent = TRUE)
      if(inherits(rval, "try-error") | is.null(rval)) rval <- nobs2(x, ...)
      return(rval)
    }
    vcov. <- function(x, ...) solve(vcov0(x) * nobs0(x))
    sandwich <- FALSE
  }

  ## score-CUSUM fluctuation process
  scus <- gefp(x, fit = NULL, order.by = order.by, vcov = vcov., scores = scores,
    decorrelate = decorrelate, sandwich = sandwich, parm = parm)

  ## set up functional if specified as character
  if(is.character(functional)) {
    functional <- tolower(functional)
    functional <- switch(functional,
      "dmax" = "dm",
      "maxlm" = "suplm",
      "mosum" = "maxmosum",
      functional
    )
    if(is.null(order.by) & functional %in% c("lmuo", "wdmo", "maxlmo")) {
      stop("'order.by' must provide a grouping of the observations")
    }
    functional <- switch(functional,
      "dm" = maxBB,
      "cvm" = meanL2BB,
      "suplm" = supLM(from = from, to = to),
      "range" = rangeBB,
      "lmuo" = catL2BB(factor(order.by)),
      "wdmo" = ordwmax(factor(order.by)),
      "maxlmo" = ordL2BB(factor(order.by), nproc = NCOL(scus$process), nobs = nobs, nrep = nrep),
      "maxmosum" = maxMOSUM(width = width),
      stop("Unknown efp functional.")
    )
  }

  ## if desired: plot test result
  if(plot) plot(scus, functional = functional, xlab = xlab, ...)
  
  ## return labeled test result
  rval <- sctest(scus, functional = functional)
  rval$data.name <- nam
  return(rval)
}

time.gefp <- function(x, ...)
{
  time(x$process, ...)
}

print.gefp <- function(x, ...)
{
  cat("\nGeneralized Empirical M-Fluctuation Process\n\n")
  cat("Call: ")
  print(x$call)
  cat("\n\n")
  cat("Fitted model: ")
  print(x$fitted.model)
  cat("\n")
  invisible(x)
}


efpFunctional <- function(functional = list(comp = function(x) max(abs(x)), time = max),
		     boundary = function(x) rep(1, length(x)),
		     computePval = NULL,
		     computeCritval = NULL,
		     plotProcess = NULL,
		     lim.process = "Brownian bridge",
		     nobs = 10000, nrep = 50000, nproc = 1:20, h = 0.5,
		     probs = c(0:84/100, 850:1000/1000))
{		     
  probs <- probs[probs != 0]

  ## store boundary function supplied
  boundary0 <- boundary

  ## compute from functional list the full functional
  ## lambda = myfun
  
  if(is.list(functional)) {
    if(identical(names(functional), c("comp", "time"))) {
      if(identical(functional[[2]], max)) {
        myfun <- function(x) {
	  rval <- apply(as.matrix(x), 1, functional[[1]])
	  functional[[2]](rval/boundary0(0:(length(rval)-1)/(length(rval)-1)))
	}
      } else
        myfun <- function(x) functional[[2]](apply(as.matrix(x), 1, functional[[1]]))
    }
    else if(identical(names(functional), c("time", "comp")))
      myfun <- function(x) functional[[2]](apply(as.matrix(x), 2, functional[[1]]))
    else  
      stop(paste(dQuote("functional"), "should be a list with elements", dQuote("comp"), "and", dQuote("time")))
  } else {
    myfun <- functional
  }

  ## setup computeStatistic function
  computeStatistic <- function(x) {
    if (all(as.matrix(x)[1,] < .Machine$double.eps)) 
      x <- as.matrix(x)[-1,]
    else x <- as.matrix(x)
    myfun(x)
  }

  ## if missing simulate values for
  ## computePval and computeCritval
  
  if(is.null(computePval) & is.null(computeCritval)) {
    if(is.null(nproc)) {
      z <- simulateBMDist(nobs = nobs, nrep = nrep, nproc = 1,
             h = h, lim.process = lim.process, functional = myfun)
      
      zquant <- c(0, quantile(z, probs = probs))
      rm(z)

      pfun <- approxfun(zquant, 1 - c(0, probs))
      computePval <- function(x, nproc = 1) {
        1 - (1 - ifelse(x > max(zquant), 0, pfun(x)))^nproc
      }
      
      cfun <- approxfun(c(0, probs), zquant)
      computeCritval <- function(alpha, nproc = 1) cfun((1-alpha)^(1/nproc))

    } else {
      z <- matrix(rep(0, nrep * length(nproc)), ncol = length(nproc))
      colnames(z) <- as.character(nproc)
      for(i in nproc)
        z[, as.character(i)] <- simulateBMDist(nobs = nobs, nrep = nrep, nproc = i,
               h = h, lim.process = lim.process, functional = myfun)

      zquant <- rbind(0, apply(z, 2, function(x) quantile(x, probs = probs)))
      rm(z)
      computePval <- function(x, nproc = 1) {
        if(as.character(nproc) %in% colnames(zquant)) {
          pfun <- approxfun(zquant[, as.character(nproc)], 1 - c(0, probs))
          ifelse(x > max(zquant[, as.character(nproc)]), 0, pfun(x))
	}
	else stop("insufficient simulated values: cannot compute p value")
      }
      computeCritval <- function(alpha, nproc = 1) {
        if(as.character(nproc) %in% colnames(zquant)) {
          cfun <- approxfun(c(0, probs), zquant[, as.character(nproc)])
          cfun(1 - alpha)
        } else stop("insufficient simulated values: cannot compute critical value")
      }
    }
  }

  if(is.null(computeCritval)) {
    nproc <- NULL
    computeCritval <- function(alpha, nproc = 1)
      uniroot(function(y) {computePval(y, nproc = nproc) - alpha}, c(0, 1000))$root
  }

  if(is.null(computePval)) {
    nproc <- NULL
    computePval <- function(x, nproc = 1)
      uniroot(function(y) {computeCritval(y, nproc = nproc) - x}, c(0, 1))$root
  }


  ## define sensible default plotting method

  if(is.null(plotProcess)) {
    if(is.list(functional)) {

      if(identical(names(functional), c("comp", "time"))) {

      ## lambda = lambda_time(lambda_comp(x))
      ## aggregate first over components then over time

        if(identical(functional[[2]], max)) {

      ## special case: lambda = max(lambda_comp(x))
      ## can also use boundary argument: b(t) = critval * boundary(t)
    
          plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
	    xlab = NULL, ylab = NULL, main = x$type.name, ylim = NULL,
	    boundary = TRUE, ...)
	  {
            n <- x$nobs
	    bound <- computeCritval(alpha = alpha, nproc = NCOL(x$process)) * boundary0(0:n/n)
	    bound <- suppressWarnings(zoo(bound, time(x)))
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
	      proc <- suppressWarnings(zoo(apply(as.matrix(x$process), 1, functional[[1]]), time(x)))
	    
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
	        if(paste(deparse(functional[[1]]), collapse = "") == "function (x) max(abs(x))") {
	          if(boundary) {
		    do.call("lines", c(list(bound), bargs))
		    do.call("lines", c(list(-bound), bargs))
		  }
	        }	      
	      }
	      plot(x$process, xlab = xlab, ylab = ylab, main = main, panel = panel, ylim = ylim, ...)
	    }
	  }

        } else {

      ## nothing specific known about lambda_time
      ## plot: first aggregate, add critval and statistic

          plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
	    xlab = NULL, ylab = NULL, main = x$type.name, ylim = NULL,
	    boundary = TRUE, statistic = TRUE, ...)
	  {
            n <- x$nobs
	    bound <- computeCritval(alpha = alpha, nproc = NCOL(x$process)) * boundary0(0:n/n)
	    bound <- suppressWarnings(zoo(bound, time(x)))
            stat <- computeStatistic(x$process)
	    stat <- suppressWarnings(zoo(rep(stat, length(time(x))), time(x)))
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
	      proc <- suppressWarnings(zoo(apply(as.matrix(x$process), 1, functional[[1]]), time(x)))
	    
	      if(is.null(ylab)) ylab <- "Empirical fluctuation process"
	      if(is.null(ylim)) ylim <- range(c(range(proc), range(bound), range(stat)))
	    
	      plot(proc, xlab = xlab, ylab = ylab, main = main, ylim = ylim, ...)
	      abline(0, 0)
	      if(boundary) do.call("lines", c(list(bound), bargs))
	      if(statistic) lines(stat, lty = 2)
	    } else {
	      panel <- function(x, ...)
	      {
                lines(x, ...)
	        abline(0, 0)
	      }
	      plot(x$process, xlab = xlab, ylab = ylab, main = main, panel = panel, ...)
	    }
	  }
        }
      }
    
      else if(identical(names(functional), c("time", "comp"))) {

      ## lambda = lambda_comp(lambda_time(x))

        plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
	    xlab = NULL, ylab = NULL, main = x$type.name, xlim = NULL, ylim = NULL,
	    boundary = TRUE, ...)
        {
          k <- NCOL(x$process)
          bound <- computeCritval(alpha = alpha, nproc = NCOL(x$process)) * boundary0(1:k/k)
	  ## for pretty printing
	  bound <- c(bound[1], bound, bound[k])
          stat <- rep(computeStatistic(x$process), k+2)
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
	    proc <- apply(as.matrix(x$process), 2, functional[[1]])
	  
	    xlabels <- colnames(x$process)
	    if(is.null(xlabels)) xlabels <- paste("Series", 1:k)
	    if(is.null(xlab)) xlab <- "Component"
	    if(is.null(ylab)) ylab <- "Statistic"
            if(is.null(ylim)) ylim <- range(c(range(proc), range(bound), range(stat), 0))
            if(is.null(xlim)) xlim <- c(0.85, k + 0.15)
	    
	    plot(1:k, proc, xlab = xlab, ylab = ylab, main = main, xlim = xlim, ylim = ylim, axes = FALSE, type = "h", ...)
	    points(1:k, proc)
	    box()
	    axis(2)
	    axis(1, at = 1:k, labels = xlabels)
	    abline(0, 0)
	    if(boundary) do.call("lines", c(list(c(0.9, 1:k, k+0.1), bound), bargs))
	    if(!identical(functional[[2]], max)) lines(c(0.9, 1:k, k+0.1), stat, lty = 2)	    
	  } else {
            if(is.null(xlab)) {
	      if(!is.null(x$order.name)) xlab <- x$order.name
	        else xlab <- "Time"
            }
	    panel <- function(x, ...)
	    {
              lines(x, ...)
              abline(0, 0)
	    }
            plot(x$process, xlab = xlab, ylab = ylab, main = main, panel = panel, ...)
	  }      
        }

      }
    
    } else {

      ## lambda = lambda(x)
      ## functional is already the full functional lambda
      ## for plotting: just plot raw process
      plotProcess <- function(x, alpha = 0.05, aggregate = FALSE,
        xlab = NULL, ylab = NULL, main = x$type.name, ...)
      {
	if(is.null(xlab)) {
	  if(!is.null(x$order.name)) xlab <- x$order.name
	    else xlab <- "Time"
	}
        if(aggregate) warning("aggregation not available")
        panel <- function(x, ...) {
          lines(x, ...)
	  abline(0, 0)
        }
        plot(x$process, xlab = xlab, ylab = ylab, main = main, panel = panel, ...)
      }  
    }
  }
  
  rval <- list(plotProcess = plotProcess,
               computeStatistic = computeStatistic,
	       computePval = computePval,
	       computeCritval = computeCritval,
	       boundary = boundary0,
	       lim.process = lim.process,
	       nobs = nobs, nrep = nrep, nproc = nproc)
	       
  class(rval) <- "efpFunctional"
  return(rval)
}

simulateBMDist <- function(nobs = 1000, nrep = 5000, nproc = 1,
                         lim.process = "Brownian bridge", 
			 h = 0.5, functional = max)
{
  lim.process <- match.arg(lim.process,
    c("Brownian motion", "Brownian motion increments", 
      "Brownian bridge", "Brownian bridge increments"))
  rval <- numeric(nrep)
  
  switch(lim.process,
  
  "Brownian motion" = {
    for(i in 1:nrep) {
      x <- matrix(rnorm(nproc * nobs), ncol = nproc)
      x <- apply(x, 2, cumsum)
      x <- rbind(0, x)/sqrt(nobs)
      rval[i] <- functional(x[-1,,drop = FALSE])
    }
  },
  
  "Brownian motion increments" = {
    nh <- floor(nobs * h)
    for(i in 1:nrep) {
      x <- matrix(rnorm(nproc * nobs), ncol = nproc)
      x <- apply(x, 2, cumsum)
      x <- rbind(0, x)/sqrt(nobs)
      x <- apply(x, 2, function(z) z[-(1:nh)] - z[1:(nobs-nh+1)])
      rval[i] <- functional(x)
    }
  },
  
  "Brownian bridge" = {
    for(i in 1:nrep) {
      x <- matrix(rnorm(nproc * nobs), ncol = nproc)
      x <- apply(x, 2, function(z) cumsum(z - mean(z)))
      x <- rbind(0, x)/sqrt(nobs)
      rval[i] <- functional(x[-1,,drop = FALSE])
    }
  },
  
  "Brownian bridge increments" = {
    nh <- floor(nobs * h)
    for(i in 1:nrep) {
      x <- matrix(rnorm(nproc * nobs), ncol = nproc)
      x <- apply(x, 2, function(z) cumsum(z - mean(z)))
      x <- rbind(0, x)/sqrt(nobs)
      x <- apply(x, 2, function(z) z[-(1:nh)] - z[1:(nobs-nh+1)])
      rval[i] <- functional(x)
    }
  })
  
  return(rval)
}

