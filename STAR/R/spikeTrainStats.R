varianceTime <- function (spikeTrain,
                          CI = c(0.95, 0.99),
                          windowSizes) {

  if (!is.spikeTrain(spikeTrain)) spikeTrain <- as.spikeTrain(spikeTrain)
  
  train.l <- length(spikeTrain)
  lastTime <- spikeTrain[train.l]
  mean.rate <- train.l/(lastTime - spikeTrain[1])
  if (!missing(windowSizes)) {
    windowSizes <- sort(windowSizes)
    if (any(windowSizes <= 0)) 
      stop(paste("Elements of", deparse(substitute(windowSizes)), 
                 "should be non-negative"))
    if (max(windowSizes) >= lastTime/5) 
      warning(paste("Some elements of", deparse(substitute(windowSizes)), 
                    "are too large and will be ignored"))
    windowSizes <- windowSizes[windowSizes < lastTime/5]
  } else {
    minSize <- 5/mean.rate
    maxSize <- lastTime/10
    windowSizes <- seq(minSize, maxSize, minSize)
    windowSizes <- windowSizes[windowSizes <= maxSize]
  }
  if (length(CI) > 2) 
    CI <- CI[1:2]
  if (any(CI >= 1 | CI <= 0)) 
    stop(paste(deparse(substitute(CI)), "components should be in (0,1)"))
  dataVariance <- sapply(windowSizes,
                         function(ws) {
                           nBins <- lastTime%/%ws
                           binSeq <- (0:nBins) * ws
                           counts <- hist(spikeTrain[spikeTrain > min(binSeq) & spikeTrain < max(binSeq)],
                                          breaks = binSeq, plot = FALSE)$counts
                           s2 <- var(counts)
                           sigma2 <- mean.rate * ws
                           ciUp <- sapply(CI, function(ci)
                                          qnorm(1 - (1 - ci)/2,sigma2, sqrt((2*sigma2^2+sigma2)/nBins)))
                           ciLow <- sapply(CI, function(ci)
                                           max(qnorm((1 - ci)/2,sigma2, sqrt((2*sigma2^2+sigma2)/nBins)), 0))
                           c(s2, sigma2, ciUp, ciLow)
                         }
                         )
  if (length(CI) == 1) { 
    resutl <- list(s2 = dataVariance[1,],
                   sigma2 = dataVariance[2,],
                   ciUp = dataVariance[3,],
                   ciLow = dataVariance[4,],
                   windowSizes = windowSizes,
                   CI = CI,
                   call = match.call()
                   )
  } else {
    result <- list(s2 = dataVariance[1, ],
                   sigma2 = dataVariance[2,],
                   ciUp = dataVariance[c(3, 4), ],
                   ciLow = dataVariance[c(5,6), ],
                   windowSizes = windowSizes,
                   CI = CI,
                   call = match.call()
                   )
  }
  class(result) <- "varianceTime"
  return(result)

}

is.varianceTime <- function(obj) {

  if(!("varianceTime" %in% class(obj))) return(FALSE)
  
  componentsNames <- c("s2",
                       "sigma2",
                       "ciUp",
                       "ciLow",
                       "call", 
                       "windowSizes",
                       "CI")
  
  if (!all(names(obj) %in% componentsNames)) return(FALSE)
  TRUE
}

plot.varianceTime <- function (x,
                               style = c("default", "Ogata"), 
                               unit = "s", xlab, ylab, main,
                               sub, xlim, ylim,
                               ...) {

  varianceTimeObj <- x
  if (!is.varianceTime(varianceTimeObj))
    stop(paste(deparse(substitute(varianceTimeObj)), "is not a varianceTime object."))
  
  if (missing(xlab)) 
    xlab <- paste("Time (", unit, ")", sep = "")
  if (missing(ylab)) 
    ylab <- "Variance"
  if (missing(main)) 
    main <- "Estimated Variance-Time Curve and Theoretical Poisson"
  if (missing(xlim)) 
    xlim <- c(0, max(varianceTimeObj$windowSizes))
  if (missing(ylim)) {
    ylim <- c(0, max(c(varianceTimeObj$s2, varianceTimeObj$ciUp)) * 1.01)
  }
  if (missing(sub)) {
    if (is.null(varianceTimeObj$call[["CI"]])) {
      CItxt <- paste(eval(formals(varianceTime)$CI) * 100,collapse = " and ")
    } else {
      CItxt <- paste(eval(varianceTimeObj$call[["CI"]]) * 100, collapse = " and ")
    }
    sub <- paste("Conf. Interval at: ", CItxt, sep = "")
  }
  X <- varianceTimeObj$windowSizes
  if (style[1] == "Ogata") {
    plot(X, varianceTimeObj$s2, type = "n", xlab = xlab, 
         ylab = ylab, main = main, xlim = xlim, ylim = ylim, 
         sub = sub, xaxs = "i", yaxs = "i")
    if (is.null(dim(varianceTimeObj$ciLow))) {
      lines(X, varianceTimeObj$ciLow, lty = 2)
      lines(X, varianceTimeObj$ciUp, lty = 2)
    } else {
      apply(varianceTimeObj$ciLow, 1, function(Y) lines(X,Y, lty = 2))
      apply(varianceTimeObj$ciUp, 1, function(Y) lines(X,Y, lty = 2))
    }
    lines(X, varianceTimeObj$sigma2)
    points(X, varianceTimeObj$s2, pch = 3, ...)
  } else {
    plot(X, varianceTimeObj$s2, type = "n", xlab = xlab, 
         ylab = ylab, main = main, xlim = xlim, ylim = ylim, 
         xaxs = "i", yaxs = "i", sub = sub)
    if (is.null(dim(varianceTimeObj$ciUp))) {
      polygon(c(X, rev(X)),
              c(varianceTimeObj$ciUp,rev(varianceTimeObj$ciLow)),
              lty = 0,
              col = "grey80")
    } else {
      polygon(c(X, rev(X)),
              c(varianceTimeObj$ciUp[2,],rev(varianceTimeObj$ciLow[2, ])),
              lty = 0, 
              col = "grey30")
      polygon(c(X, rev(X)),
              c(varianceTimeObj$ciUp[1,],rev(varianceTimeObj$ciLow[1, ])),
              lty = 0, 
              col = "grey80")
    }
    lines(X, varianceTimeObj$sigma2, lwd = 2)
    lines(X, varianceTimeObj$s2, col = 2, lwd = 2, ...)
  }
  
}


acf.spikeTrain <- function (spikeTrain,
                            lag.max = NULL,
                            type = c("correlation", "covariance", "partial"),
                            plot = TRUE,
                            na.action = na.fail, 
                            demean = TRUE,
                            xlab = "Lag (in isi #)",
                            ylab = "ISI acf", 
                            main,
                            ...) {

  if (!is.spikeTrain(spikeTrain)) spikeTrain <- as.spikeTrain(spikeTrain)
  isi <- diff(spikeTrain)
  if (missing(main)) 
    main <- paste("Train", deparse(substitute(spikeTrain)), "ISI acf")
  acf(isi, lag.max = lag.max, type = type, plot = plot, na.action = na.action, 
      demean = demean, xlab = xlab, main = main, ylab = ylab, 
      ...)
  
}

renewalTestPlot <- function (spikeTrain,
                             lag.max=NULL,
                             d=max(c(2,sqrt(length(spikeTrain)) %/% 5)),
                             orderPlotPch=ifelse(length(spikeTrain)<=600,1,"."), 
                             ...) {

  if (!is.spikeTrain(spikeTrain)) spikeTrain <- as.spikeTrain(spikeTrain)
  if (length(spikeTrain) < 50) 
    stop(paste(deparse(substitute(spikeTrain)), "contains less than 50 events."))
  m <- matrix(c(1:4), nrow = 2, ncol = 2, byrow = TRUE)
  oldpar <- par(mar = c(5, 4, 2, 2))
  layout(m)
  on.exit(par(oldpar))
  isi <- diff(spikeTrain)
  isi.o <- rank(isi)/length(isi)
  isi.l <- length(isi)
  if (is.null(lag.max)) 
    lag.max <- round(10 * log10(isi.l))
  lag.max <- min(isi.l - 1, lag.max)
  grid <- seq(0, 1, length.out = d + 1)
  getChi2 <- function(lag) {
    isi.x <- isi.o[1:(isi.l - lag.max)]
    isi.y <- isi.o[(1 + lag):(isi.l - lag.max + lag)]
    isi.x <- as.integer(cut(isi.x, breaks = grid))
    isi.y <- as.integer(cut(isi.y, breaks = grid))
    counts <- matrix(0, nrow = d, ncol = d)
    for (i in seq(along.with = isi.x))
      counts[isi.x[i], isi.y[i]] <- counts[isi.x[i], isi.y[i]] + 1
    chisq.test(counts, ...)
  }
  chi2seq <- lapply(1:lag.max, getChi2)
  minChi2 <- qchisq(0.025, df = chi2seq[[1]]$parameter)
  maxChi2 <- qchisq(0.975, df = chi2seq[[1]]$parameter)
  chi2V <- sapply(chi2seq, function(l) l$statistic)
  outOf95 <- chi2V < minChi2 | chi2V > maxChi2

  plot(isi.o[-length(isi.o)] * isi.l, isi.o[-1] * isi.l, xlab = quote(O[k]), 
       ylab = quote(O[k + 1]), main = "Order Statistic Correlation at Lag 1", 
       type = "n")

  sapply(grid[-c(1, d + 1)],
         function(idx) {
           abline(v = idx * isi.l, col = "grey")
           abline(h = idx * isi.l, col = "grey")
         }
         )

  points(isi.o[-length(isi.o)] * isi.l, isi.o[-1] * isi.l, 
         pch = orderPlotPch)
  plot(isi.o[-(0:-1 + length(isi.o))] * isi.l,
       isi.o[-(1:2)] * isi.l,
       xlab = quote(O[k]), ylab = quote(O[k + 2]),
       main = "Order Statistic Correlation at Lag 2", 
       type = "n")
  
  sapply(grid[-c(1, d + 1)],
         function(idx) {
           abline(v = idx * isi.l, col = "grey")
           abline(h = idx * isi.l, col = "grey")
         }
         )
  
  points(isi.o[-(0:-1 + length(isi.o))] * isi.l,
         isi.o[-(1:2)] * isi.l,
         pch = orderPlotPch)
  
  acf.spikeTrain(spikeTrain,
                 lag.max = lag.max, main = "")

  plot(1:lag.max, chi2V, type = "n",
       xlab = "Lag (in isi #)", 
       ylab = quote(chi^2),
       main = expression(paste(chi^2, " Statistics")), 
       xlim = c(0, lag.max),
       ylim = c(min(qchisq(0.01, df = chi2seq[[1]]$parameter), min(chi2V)),
         max(qchisq(0.99, df = d^2 - 2 * (d - 1) - 1), max(chi2V)))
       )
  
  polygon(c(0, 0, lag.max + 1, lag.max + 1),
          c(minChi2, maxChi2, maxChi2, minChi2),
          col = "grey80", lty = 0
          )

  points(1:lag.max, chi2V,
         pch = c(16, 17)[outOf95 + 1],
         col = c("black", "grey30")[outOf95 + 1]
         )

}
