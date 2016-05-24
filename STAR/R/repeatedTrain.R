as.repeatedTrain <- function(x) {
  if (!is.list(x)) 
    x <- list(x)
  if (!is.null(names(x))) {
    xN <- names(x)
  } else {
    xN <- paste("stim.",1:length(x))
  }
  
  x <- lapply(x, as.spikeTrain)
  names(x) <- xN
  class(x) <- "repeatedTrain"
  x
}

is.repeatedTrain <- function(x) {
  if (!("repeatedTrain" %in% class(x))) return(FALSE)
  if (!all(sapply(x,is.spikeTrain))) return(FALSE)
  return(TRUE)
}

raster <- function(x,
                   stimTimeCourse = NULL,
                   colStim = "grey80", 
                   xlim, pch, xlab, ylab, main,
                   ...) {

  do.call(plot.repeatedTrain,as.list(match.call()[-1]))

}

plot.repeatedTrain <- function (x,
                                stimTimeCourse = NULL,
                                colStim = "grey80", 
                                xlim, pch, xlab, ylab, main,
                                ...) {
  if (!is.repeatedTrain(x)) x <- as.repeatedTrain(x)
  nbTrains <- length(x)
  if (missing(xlim)) 
    xlim <- c(0, ceiling(max(sapply(x, max))))
  if (missing(xlab)) 
    xlab <- "Time (s)"
  if (missing(ylab)) 
    ylab <- "trial"
  if (missing(main)) 
    main <- paste(deparse(substitute(x)), "raster")
  if (missing(pch)) 
    pch <- ifelse(nbTrains <= 20, "|", ".")
  oldpar <- par(mar = c(5, 4, 2, 1))
  on.exit(par(oldpar))
  acquisitionDuration <- max(xlim)
  plot(c(0, acquisitionDuration), c(0, nbTrains + 1), type = "n", 
       xlab = xlab, ylab = ylab, xlim = xlim, ylim = c(1, nbTrains + 
                                                1), bty = "n", main = main, ...)
  if (!is.null(stimTimeCourse)) {
    rect(stimTimeCourse[1], 0.1, stimTimeCourse[2], nbTrains + 
         0.9, col = colStim, lty = 0)
  }
  invisible(sapply(1:nbTrains, function(idx) points(x[[idx]], 
                                                    numeric(length(x[[idx]])) + idx, pch = pch)))
  axis(2, at = 1:nbTrains)
  axis(1)

}

plot.psth <- function(x,
                      stimTimeCourse = NULL,
                      colStim = "grey80", 
                      colCI = NULL,
                      xlab, ylab, main, xlim, 
                      ylim, lwd = 2, col = 1,
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

  withCI <- !is.null(x$ciUp)

  originalBreaksArg <- x$call[["breaks"]]
  if (!is.null(originalBreaksArg) && (length(eval(originalBreaksArg))==2))
    smoothOne <- TRUE
  else
    smoothOne <- FALSE
  if (missing(xlim)) {
    if (!smoothOne) xlim <- range(x$breaks)
    else xlim <- range(x$mids) + c(-0.5,0.5)*x$breaks["step"]
  }  
  if (missing(ylim))
    ylim <- c(0, ifelse(withCI, max(x$ciUp), max(x$freq))) ## Bug fix thanks to Gregory Jefferis

  

  plot(xlim,
       ylim,
       type = "n", 
       xlab = xlab,
       ylab = ylab,
       main = main,
       xlim = xlim, 
       ylim = ylim,
       ...)
  
  if (!is.null(stimTimeCourse)) {
    rect(stimTimeCourse[1], 0, stimTimeCourse[2], ylim[2], 
         col = colStim, lty = 0)
  }
  if (withCI) {
    if (is.null(colCI)) {
      if (!smoothOne) {
        sapply(1:length(x$mids),
               function(idx)
               segments(x$breaks[idx], x$ciLow[idx],
                        x$breaks[idx + 1], x$ciLow[idx], 
                        lty = 2)
               )
        sapply(1:length(x$mids),
               function(idx)
               segments(x$breaks[idx], x$ciUp[idx],
                        x$breaks[idx + 1], x$ciUp[idx], 
                        lty = 2)
               )
      } else {
        lines(x$mids, x$ciLow, lty = 2)
        lines(x$mids, x$ciUp, lty = 2)
      }
    } else {
      if (!smoothOne) {
        sapply(1:length(x$mids),
               function(idx)
               rect(x$breaks[idx], x$ciLow[idx],
                    x$breaks[idx + 1], x$ciUp[idx], 
                    lty = 0, col = colCI))
      } else {
        polygon(c(x$mids, rev(x$mids)), c(x$ciLow, rev(x$ciUp)), 
                col = colCI, border = NA)
      }
    }
  }
  if (!smoothOne) {
    invisible(sapply(1:length(x$mids),
                     function(idx) {
                       segments(x$breaks[idx], x$freq[idx],
                                x$breaks[idx + 1], x$freq[idx],
                                lwd = lwd, col = col)
                       if (idx == 1)
                         segments(x$breaks[idx], 0, x$breaks[idx], 
                                  x$freq[idx], lwd = lwd, col = col)
                       if (idx != length(x$mids)) {
                         segments(x$breaks[idx + 1],
                                  x$freq[idx],
                                  x$breaks[idx + 1],
                                  x$freq[idx + 1],
                                  lwd = lwd,
                                  col = col)
                       } else {
                         segments(x$breaks[idx + 1],
                                  x$freq[idx],
                                  x$breaks[idx + 1],
                                  0, lwd = lwd, col = col)
                       }
                     }
                     )
              )
  } else {
    lines(x$mids, x$freq, col = col, lwd = lwd)
  }
  abline(h = 0)
}

psth <- function (repeatedTrain,
                  breaks=20,
                  include.lowest = TRUE,
                  right = TRUE, 
                  plot = TRUE,
                  CI = 0.95,
                  ...
                  ) {
  
  if (!is.repeatedTrain(repeatedTrain))
    repeatedTrain <- as.repeatedTrain(repeatedTrain)

  if (!inherits(breaks,"numeric"))
    stop("breaks should be a numeric.")
 
  nbTrials <- length(repeatedTrain)
  breaksName <- deparse(substitute(breaks))
  l <- floor(min(sapply(repeatedTrain,function(l) l[1])))
  r <- ceiling(max(sapply(repeatedTrain,function(l) l[length(l)])))
      
  if (length(breaks) != 2) {
    if (length(breaks)==1) breaks <- seq(l,r,length.out=breaks+1)
    counts <- t(sapply(repeatedTrain,
                       function(train)
                       hist(x = unclass(train), breaks = breaks,
                            include.lowest = include.lowest, 
                            right = right, plot = FALSE)$counts
                       )
                )
    h <- list(breaks = breaks,
              counts = counts,
              mids = breaks[-length(breaks)]+diff(breaks)/2)
    f <- colSums(counts)/(nbTrials * diff(h$breaks))
  } else {
    if (is.null(names(breaks))) {
      bw <- breaks[1]
      step <- breaks[2]
    } else {
      if (!all(names(breaks) %in% c("bw", "step"))) 
        stop(paste(breaksName, "should have named elements: bw and step"))
      bw <- as.numeric(breaks["bw"])
      step <- as.numeric(breaks["step"])
    } ## End of conditional on is.null(names(breaks)
    bwh <- bw/2
    breaks <- c(bw = bw, step = step)
    mids <- seq(bwh, r - bwh, by = step)
    counts <- t(sapply(repeatedTrain,
                       function(train)
                       sapply(mids,
                              function(m)
                              ifelse(right,
                                     sum(m - bwh < train & train <= m + bwh),
                                     sum(m - bwh <= train & train < m + bwh))
                              )
                       )
                )
    h <- list(breaks = breaks,
              counts = counts,
              mids = mids)
    f <- colSums(counts)/(nbTrials * bw)
  } ## End of conditional on length(breaks) != 2
  if (!is.null(CI)) {
    expectedCount <- colSums(counts)
    ciUp <- sapply(1:length(h$mids),
                   function(idx) qpois(1-(1-CI)/2,expectedCount[idx])
                   )
    if (length(breaks) > 2) ciUp <- ciUp / (nbTrials * diff(h$breaks))
    else ciUp <- ciUp / (nbTrials * bw)
    ciLow <- sapply(1:length(h$mids),
                    function(idx) qpois((1-CI)/2,expectedCount[idx])
                    )
    if (length(breaks) > 2) ciLow <- ciLow / (nbTrials * diff(h$breaks))
    else ciLow <- ciLow / (nbTrials * bw)
  }

  if (!is.null(CI)) {
    result <- list(freq = f, ciUp = ciUp, ciLow = ciLow, 
                   breaks = h$breaks, mids = h$mids, counts = h$counts, 
                   nbTrials = nbTrials, call = match.call())
    class(result) <- "psth"
  } else {
    result <- list(freq = f, ciUp = NULL, ciLow = NULL, 
                   breaks = h$breaks, mids = h$mids, counts = h$counts, 
                   nbTrials = nbTrials, call = match.call())
    class(result) <- "psth"
  }
  
  if (plot) {
    plot(result,...)
  } else {
    return(result)
  }
  
}

print.repeatedTrain <- function(x,...) {
  main <- "Raster plot"
  plot(x,main=main,...)
}

summary.repeatedTrain <- function(object,
                                  responseWindow,
                                  acquisitionWindow,
                                  ...) {
  repeatedTrain <- object
  rm(object)
  nbRepeates <- length(repeatedTrain)
  if (missing(acquisitionWindow))
    acquisitionWindow <- c(min(sapply(repeatedTrain,function(l) floor(l[1]))),
                           max(sapply(repeatedTrain,function(l) ceiling(l[length(l)])))
                           )
  
  acquisitionDuration <- diff(acquisitionWindow)
  if (!missing(responseWindow)) {
    if (length(responseWindow) != 2) stop("Wrong responseWindow.")
    if ((responseWindow[1] < acquisitionWindow[1]) ||
        (responseWindow[2] > acquisitionWindow[2]) )
      warning("responseWindow does not make much sense.")
    stimDuration <- diff(responseWindow)
  } else {
    responseWindow <- NULL
    stimDuration <- NULL
  }
  theStats <- sapply(repeatedTrain,
                     function(l) {
                       nb <- length(l)
                       nu <- nb/acquisitionDuration
                       result <- c(nb=nb,nu=nu)
                       if (!is.null(stimDuration)) {
                         goodOnes <- responseWindow[1] < l & l < responseWindow[2]
                         nbR <- sum(goodOnes)
                         nuR <- nbR/stimDuration
                         result <- c(result,c(nbR=nbR,nuR=nuR))
                       }
                       result
                     }
                     )
  stationarityTest <- chisq.test(theStats["nb",],...)$p.value
  if (!is.null(stimDuration)) {
    stationarityTestR <- chisq.test(theStats["nbR",],...)$p.value
  } else {
    stationarityTestR <- NULL
  }

  result <- list(nbRepeates=nbRepeates,
                 acquisitionWindow=acquisitionWindow,
                 responseWindow=responseWindow,
                 stats=theStats,
                 globalPval=stationarityTest,
                 responsePval=stationarityTestR)
  class(result) <- "summary.repeatedTrain"
  return(result)
}

print.summary.repeatedTrain <- function(x,...) {

  obj <- x
  rm(x)
  cat(paste(obj$nbRepeates, " repeats in the object.\n"))
  print(t(obj$stats))
  cat(paste("The p value of the chi 2 test for the stationarity accross repeats is:\n",
            obj$globalPval,
            ".\n",sep="")
      )
  if (!is.null(obj$responseWindow))
    cat(paste("The p value of the chi 2 test for the stationarity accross repeats during the stim. is:\n",
            obj$responsePval,
            ".\n",sep="")
      )
  
}


df4counts <- function(repeatedTrain,
                      breaks=length(repeatedTrain)
                      ) {
  
  if (!is.repeatedTrain(repeatedTrain))
    repeatedTrain <- as.repeatedTrain(repeatedTrain)

  nbTrials <- length(repeatedTrain)
  breaksName <- deparse(substitute(breaks))
  trainNames <- names(repeatedTrain)
  if (is.null(trainNames))
    trainNames <- paste("trial",1:nbTrials)
  
  l <- floor(min(sapply(repeatedTrain,function(l) l[1])))
  r <- ceiling(max(sapply(repeatedTrain,function(l) l[length(l)])))

  if (length(breaks) == 1) {
    ## breaks is interpreted as a number of bins
    breaks <- seq(l,r,length.out=breaks+1)
  } ## End of conditional on length(breaks) == 1

  counts <- sapply(repeatedTrain,
                   function(train)
                   hist(x=train[l <= train & train <= r],
                        breaks=breaks, plot=FALSE)$counts
                   )

  mids <- breaks[-length(breaks)]+diff(breaks)/2
  nb <- length(mids)
  data.frame(Count=as.numeric(counts),
             Bin=factor(rep(1:nb,nbTrials)),
             Trial=factor(rep(trainNames,each=nb)),
             Rate=as.numeric(counts)/diff(breaks),
             Time=rep(mids,nbTrials)
             )

}
