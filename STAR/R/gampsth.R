gampsth <- function (repeatedTrain,
                     binSize=0.025,
                     k=100,
                     bs="tp",
                     plot=TRUE,
                     ...) {
  
  if (!is.repeatedTrain(repeatedTrain)) 
    repeatedTrain <- as.repeatedTrain(repeatedTrain)

  bigTrain <- sort(unlist(repeatedTrain))
  minSpikeTime <- floor(min(bigTrain))
  maxSpikeTime <- ceiling(max(bigTrain))
  Time <- seq(minSpikeTime,maxSpikeTime,by=binSize)
  if (max(bigTrain) > max(Time)) Time <- c(Time,max(Time)+binSize)
  H <- hist(bigTrain,breaks=Time,plot=FALSE)
  Count <- H$counts
  Time <- H$mids

  PoissonF <- gam(Count ~ s(Time,k=k,bs=bs),family=poisson)
  nbTrials <- length(repeatedTrain)

  Y <- predict(PoissonF,type="response",se.fit=TRUE)
  Y$fit <- Y$fit/nbTrials/binSize
  Y$se.fit <- Y$se.fit/nbTrials/binSize

  ## create a function returning the intensity
  lambdaFct <- function(t) {
    as.numeric(predict(PoissonF,
                       newdata=data.frame(Time=t),
                       type="response")/nbTrials/binSize
               )
  }

  ## create a function returning the integrated intensity
  LambdaFct <- function(t) {
    t <- sort(t)
    t <- t[minSpikeTime <= t & t <= maxSpikeTime]
    tL <- length(t)
    result <- numeric(tL)
    result[1] <- integrate(lambdaFct,minSpikeTime,t[1])$value
    if (tL > 1) for (i in 2:tL) result[i] <- integrate(lambdaFct,t[i-1],t[i])$value
    cumsum(result)
  }

  result <- list(freq=Y$fit,
                 ciUp=Y$fit+1.96*Y$se.fit,
                 ciLow=Y$fit-1.96*Y$se.fit,
                 breaks=c(minSpikeTime,maxSpikeTime), 
                 mids=Time,
                 counts=Count,
                 nbTrials=nbTrials,
                 lambdaFct=lambdaFct,
                 LambdaFct=LambdaFct,
                 call=match.call())
  class(result) <- "gampsth"
    
  if (plot) {
    plot(result, ...)
  }
  else {
    return(result)
  }
  
}


summary.gampsth <- function(object, ...) {
  summary(gamObj(object))
}

print.gampsth <- function(x, ...) {
  print(gamObj(x))
}


plot.gampsth <- function (x,
                        stimTimeCourse=NULL,
                        colStim="grey80",
                        colCI=NULL, 
                        xlab, ylab, main, xlim, ylim,
                        lwd=2, col=1,
                        ...) {
  
    if (!is.null(stimTimeCourse)) {
      if (length(stimTimeCourse) != 2) 
        stop(paste(deparse(substitute(stimTimeCourse)),
                   "should be a vector with 2 elements.")
             )
    }
    if (missing(xlab)) 
      xlab <- "Time (s)"
    if (missing(ylab)) 
      ylab <- "Freq (Hz)"
    if (missing(main)) {
      nameList <- deparse(x$call[["repeatedTrain"]])
      main <- paste(nameList, "PSTH")
    }
    
    withCI <- TRUE
    smoothOne <- TRUE
    if (missing(xlim)) xlim <- x$breaks
    if (missing(ylim)) ylim <- c(0, max(x$ciUp))

    ## Generate the frame of the plot with "decorations"
    plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab, main = main, 
         xlim = xlim, ylim = ylim, ...)

    ## if a stim time course is given make it appear
    if (!is.null(stimTimeCourse)) {
      rect(stimTimeCourse[1], 0, stimTimeCourse[2], ylim[2], 
           col = colStim, lty = 0)
    }

    ## Add confidence bands
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
