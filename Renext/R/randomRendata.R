##*************************************************************************
## Censored Poisson distribution.
##
##*************************************************************************

.rpoisCens <- function(n, lambda = 1, r = 0) {
  
  P <- ppois(r, lambda = lambda, lower.tail = FALSE)
  sim <- rep(NA, n)
  
  for (i in 1L:n) {
    U <- runif(1)
    k <- r
    Pk <- 0
    while (Pk / P < U) {
      pk <- dpois(k, lambda = lambda)
      Pk <- Pk + pk
      k <- k + 1
    }
    sim[i] <- k-1
  }
  sim
}

##*************************************************************************
## Generate a 'Rendata' object.
##
##*************************************************************************

rRendata <- function(lambda = 1,
                     threshold = 0,
                     effDuration = 100,
                     distname.y = "exp",
                     par.y = c(rate = 1),
                     start = "1913-01-01",
                     name = NULL,
                     varName = "X",
                     varUnit = "?",
                     simDate = TRUE,
                     roundDate = FALSE,
                     MAX.effDuration = NULL,
                     MAX.r = rep(1L, length(MAX.effDuration)),
                     OTS.effDuration = NULL,
                     OTS.threshold = NULL) {
  
  if (distname.y == "exponential") distname.y <- "exp"
  if (is.null(name)) {
    name <- sprintf("simulated \"%s\" data", distname.y)
  }
  
  start <- as.POSIXct(start)

  res <- list(info = list(name = name,
                shortLab = name, longLab = name,
                varName = varName, varShortLab = varName, varUnit = varUnit,
                describe = sprintf("threshold = %7.3f, distname.y = \"%s\"",
                  threshold, distname.y)))
  
  ##=======================================================================
  ## if the effective
  ##=======================================================================

  if (effDuration > 0) {
    
    n <- rpois(n = 1, lambda = effDuration * lambda)
    T <- sort(runif(n, min = 0, max = effDuration))
    T <- start + T * 86400 * 365
    end <- start + effDuration * 86400 * 365
    
    nm <- paste("r", distname.y, sep = "")
    args <- c(n = n,  as.list(par.y))
    y <- do.call(nm, args)
    
    if (!is.na(threshold) && is.finite(threshold)) {
      x <- threshold + y
    } else x <- y
    
    OTdata <- data.frame(T, x)
    colnames(OTdata) <- c("date", varName)
    
    res$OTinfo <- list(start = start, end = end, effDuration = effDuration,
                       threshold = threshold)
    res$OTdata <-  OTdata
  } else {
    res$OTinfo <- NULL
    res$OTdata <- NULL
  }

  ##=======================================================================
  ## MAX part, if needed
  ##=======================================================================
  if (!is.null(MAX.effDuration)) {
    if (any(is.na(MAX.effDuration) | MAX.effDuration <= 0))
      stop("'MAX.effDuration' must contain positive numeric values")

    if (length(MAX.r) != length(MAX.effDuration))
        stop("'MAX.r' must have the same length as 'MAX.effDuration'")

    if (any(MAX.r < 1)) stop("'MAX.r' must contain integers >= 1")
    
    MAX.tot <- sum(MAX.effDuration)
    nBlocks <- length(MAX.effDuration)
    
    nm <- paste("r", distname.y, sep = "")
    breaks <- -MAX.tot + c(0, cumsum(MAX.effDuration))
    ## simDate <- FALSE

    blocks <- integer(0)
    Ts <- numeric(0)
    ys <- numeric(0)
    
    for (i in 1L:nBlocks) {
      ri <-  as.integer(MAX.r[i])
      ni <- .rpoisCens(n = 1, lambda = lambda * MAX.effDuration[i], r = ri)
      args <- c(n = ni,  as.list(par.y))
      yi <- do.call(nm, args)
      so <- sort(yi, decreasing = TRUE, index.return = TRUE)
      yi <- so$x[1L:ri]
      
      if (simDate) {
        Ti <- sort(runif(ni, min = breaks[i], max = breaks[i+1]))
        Ti <- Ti[so$ix[1L:ri]]
      } else {
        Ti <- rep(NA, ri)
      }
      Ts <- c(Ts, Ti) 
      blocks <- c(blocks, rep(i, ri))
      ys <- c(ys, yi)
    }

    if (simDate) Ts <- start + Ts * 365 * 86400

    ## build MAXdata
    MAXdata <- data.frame(block = blocks, date = Ts,
                          y = ys, comment = rep ("", length(blocks)))
    colnames(MAXdata) <- c("block", "date", varName, "comment")

    ## build MAXinfo
    start.MAX <- start + breaks * 365 * 86400
    end.MAX <- start.MAX[2:(1L + nBlocks)]
    start.MAX <- start.MAX[1L:nBlocks]
    res$MAXinfo <- data.frame(start = start.MAX,
                              end = end.MAX,
                              duration = MAX.effDuration)
    res$MAXdata <- MAXdata
    start <- start - MAX.tot * 365 * 86400
                          
  }

  ##=======================================================================
  ## OTS part, if needed
  ##=======================================================================
  
  if (!is.null(OTS.effDuration)) {
    if (any(is.na(OTS.effDuration) | OTS.effDuration <= 0))
      stop("'OTS.effDuration' must contain positive numeric values")
    if (length(OTS.threshold) != length(OTS.effDuration))
      stop("'OTS.threshold' must have the same length as 'OTS.effDuration'")

    if (any(OTS.threshold < threshold))
      stop("all elements of 'OTS.threshold' must be >= the value of 'threshold'")
    
    OTS.tot <- sum(OTS.effDuration)
    nBlocks <- length(OTS.effDuration)
    
    nm <- paste("r", distname.y, sep = "")
    breaks <- -OTS.tot + c(0, cumsum(OTS.effDuration))
    ## simDate <- FALSE
    
    blocks <- integer(0)
    Ts <- numeric(0)
    ys <- numeric(0)
    OTS.r <- rep(NA, nBlocks)
    
    for (i in 1L:nBlocks) {
      ni <- rpois(n = 1, lambda = lambda * OTS.effDuration[i])
      args <- c(n = ni,  as.list(par.y))
      yi <- do.call(nm, args)
      indi <- (yi >= OTS.threshold[i])
      yi <- yi[indi]
      OTS.r[i] <- sum(indi)

      if (OTS.r[i] > 0 ) {
        if (simDate) {
          Ti <- sort(runif(ni, min = breaks[i], max = breaks[i+1]))
          Ti <- Ti[indi]
        } else {
          Ti <- rep(NA, OTS.r[i])
        }
        Ts <- c(Ts, Ti) 
        blocks <- c(blocks, rep(i, OTS.r[i]))
        ys <- c(ys, yi)
      }
    }
    
    if (simDate) Ts <- start + Ts * 365 * 86400

    ## build OTSdata
    OTSdata <- data.frame(block = blocks, date = Ts,
                          y = ys, comment = rep("", length(blocks)))
    colnames(OTSdata) <- c("block", "date", varName, "comment")

    ## build OTSinfo
    start.OTS <- start + breaks * 365 * 86400
    end.OTS <- start.OTS[2:(1L + nBlocks)]
    start.OTS <- start.OTS[1L:nBlocks]
    res$OTSinfo <- data.frame(start = start.OTS,
                              end = end.OTS,
                              duration = OTS.effDuration,
                              threshold = OTS.threshold,
                              r = OTS.r)
    res$OTSdata <- OTSdata

  }

  class(res) <- "Rendata"
  res
  
  
  
}
                           
