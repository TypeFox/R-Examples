.mkDiscreteTrain <- function(repeatedTrain,
                             binSize
                             ) {
  
  if (!is.repeatedTrain(repeatedTrain)) 
    repeatedTrain <- as.repeatedTrain(repeatedTrain)
  bigTrain <- sort(unlist(repeatedTrain))
  minSpikeTime <- floor(min(bigTrain))
  maxSpikeTime <- ceiling(max(bigTrain))
  Time <- seq(minSpikeTime, maxSpikeTime, by = binSize)
  if (max(bigTrain) > max(Time)) 
    Time <- c(Time, max(Time) + binSize)
  H <- hist(bigTrain, breaks = Time, plot = FALSE)
  result <- data.frame(Count=H$counts,
                       Time=H$mids)
  attr(result,"minSpikeTime") <- minSpikeTime
  attr(result,"maxSpikeTime") <- maxSpikeTime
  return(result)
}

.postFitFormat <- function(gfit,
                           repeatedTrain,
                           binSize,
                           minSpikeTime,
                           maxSpikeTime
                           ) {

  nbTrials <- length(repeatedTrain)
  Time <- gfit$mf$Time
  Count <- gfit$mf$Count
  Y <- predict(gfit,
               newdata=data.frame(Time=Time),
               se.fit = TRUE)
  lambdaFct <- function(t) {
    exp(as.numeric(predict(gfit, newdata = data.frame(Time = t))))/nbTrials/binSize
  }
  LambdaFct <- function(t) {
    t <- sort(t)
    t <- t[minSpikeTime <= t & t <= maxSpikeTime]
    tL <- length(t)
    result <- numeric(tL)
    result[1] <- integrate(lambdaFct, minSpikeTime, t[1])$value
    if (tL > 1) 
      for (i in 2:tL) result[i] <- integrate(lambdaFct, 
                                             t[i - 1], t[i])$value
    cumsum(result)
  }
  result <- list(freq = exp(Y$fit)/nbTrials/binSize,
                 ciUp = exp(Y$fit + 1.96 * Y$se.fit)/nbTrials/binSize, 
                 ciLow = exp(Y$fit - 1.96 * Y$se.fit)/nbTrials/binSize,
                 breaks = c(minSpikeTime, 
                   maxSpikeTime),
                 mids = Time,
                 counts = Count,
                 nbTrials = nbTrials, 
                 lambdaFct = lambdaFct,
                 LambdaFct = LambdaFct,
                 call = match.call())
  return(result)
  
}

gsspsth0 <- function(repeatedTrain,
                     binSize = 0.025, 
                     plot = FALSE,
                     ...) {
  
  df <- .mkDiscreteTrain(repeatedTrain,
                         binSize)
  PoissonF <- gssanova0(Count ~ Time, family = "poisson", data=df, ...)
  result <- .postFitFormat(PoissonF,repeatedTrain,binSize,
                           attr(df,"minSpikeTime"),
                           attr(df,"maxSpikeTime"))
  class(result) <- "gsspsth0"
  if (plot) {
    plot(result, colCI=2)
  }
  else {
    return(result)
  }
}


gsspsth <- function(repeatedTrain,
                    binSize = 0.025, 
                    plot = FALSE,
                    ...) {
  
  df <- .mkDiscreteTrain(repeatedTrain,
                         binSize)
  PoissonF <- gssanova(Count ~ Time, family = "poisson", data = df,
                       nbasis=dim(df)[1],...)
  result <- .postFitFormat(PoissonF,repeatedTrain,binSize,
                           attr(df,"minSpikeTime"),
                           attr(df,"maxSpikeTime"))
  class(result) <- "gsspsth"
  if (plot) {
    plot(result, colCI=2)
  }
  else {
    return(result)
  }
}

.plot_spsth <- function(x,
                        stimTimeCourse,
                        colStim,
                        colCI, 
                        xlab, ylab, main, xlim, ylim,
                        lwd, col,
                        ...) {
  
  plot(xlim, ylim, type = "n",
       xlab = xlab, ylab = ylab, main = main, 
       xlim = xlim, ylim = ylim,
       ...)
  
  if (!is.null(stimTimeCourse)) {
    rect(stimTimeCourse[1], 0, stimTimeCourse[2], ylim[2], 
         col = colStim, lty = 0)
  }
  if (is.null(colCI)) {
    lines(x$mids, x$ciLow, lty = 2)
    lines(x$mids, x$ciUp, lty = 2)
  } else {
    polygon(c(x$mids, rev(x$mids)), c(x$ciLow, rev(x$ciUp)), 
            col = colCI, border = NA)
  }
  lines(x$mids, x$freq, col = col, lwd = lwd)
  abline(h = 0)
  
}

plot.gsspsth <- function(x,
                         stimTimeCourse = NULL,
                         colStim = "grey80",
                         colCI = NULL, 
                         xlab, ylab, main, xlim, ylim,
                         lwd = 2, col = 1,
                         ...) {
  
  if (!is.null(stimTimeCourse)) {
    if (length(stimTimeCourse) != 2) 
      stop(paste(deparse(substitute(stimTimeCourse)), "should be a vector with 2 elements."))
  }
  if (missing(xlab)) xlab <- "Time (s)"
  if (missing(ylab)) ylab <- "Freq (Hz)"
  if (missing(main)) {
    nameList <- deparse(x$call[["repeatedTrain"]])
    main <- paste(nameList, "PSTH")
  }
  if (missing(xlim)) 
    xlim <- x$breaks
  if (missing(ylim)) 
    ylim <- c(0, max(x$ciUp[!is.na(x$ciUp)]))

  .plot_spsth(x,stimTimeCourse,colStim,colCI, 
              xlab, ylab, main, xlim, ylim, lwd, col,
              ...)
}

plot.gsspsth0 <- plot.gsspsth


gssObj <- function(object,
                   ...) {

  UseMethod("gssObj")
}

gssObj.gsspsth <- function(object,...) {

  evalq(gfit, envir=environment(object$lambdaFct))

}


summary.gsspsth <- function(object, ...) {
  summary(gssObj(object))
}

print.gsspsth <- function(x, ...) {
  print(gssObj(x))
}

gssObj.gsspsth0 <- gssObj.gsspsth
summary.gsspsth0 <- summary.gsspsth
print.gsspsth0 <- print.gsspsth


simulate.gsspsth0 <- function(object,
                              nsim=1,
                              seed=NULL,
                              ...) {
  
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)                       ## initialize the RNG if necessary
  if(is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  ## find out the largest frequency
  idx.max <- which.max(object$freq)
  if (idx.max == 1) {
    lwr <- 1
    upr <- 2
  } else {
    freq.length <- length(object$freq)
    if (idx.max == freq.length) {
      lwr <- freq.length-1
      upr <- freq.length
    } else {
      lwr <- idx.max-1
      upr <- idx.max+1
    }
  }
  lwr <- object$mids[lwr]
  upr <- object$mids[upr]
  f.max <- optimize(object$lambdaFct,lower=lwr,upper=upr,maximum=TRUE)$objective

  time.range <- object$breaks
  
  mk.train <- function(f.max,time.range) {
    first.guess <- diff(time.range)*f.max*2
    result <- cumsum(rexp(first.guess,f.max))+time.range[1]
    result.max <- result[first.guess]
    while (result.max < time.range[2]) {
      result <- c(result,cumsum(rexp(first.guess,f.max))+result.max)
      first.guess <- length(result)
      result.max <- result[first.guess]
    }
    result <- result[result < time.range[2]]
    ratio <- object$lambdaFct(result)/f.max
    result <- result[runif(length(result))<=ratio]
    class(result) <- "spikeTrain"
    result
  }
  
  result <- lapply(1:nsim,
                   function(sIdx) {
                     ans <- lapply(1:object$nbTrials,
                                   function(tIdx) mk.train(f.max,time.range)
                                   )
                     class(ans) <- "repeatedTrain"
                     ans
                   }
                   )
  attr(result, "seed") <- RNGstate

  if (nsim==1) result <- result[[1]]

  result
}

simulate.gsspsth <- simulate.gsspsth0
