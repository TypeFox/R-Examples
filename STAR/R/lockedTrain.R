lockedTrain <- function(stRef,
                        stTest,
                        laglim,
                        acquisitionWindow
                        ) {

  if (missing(stTest)) {
    stTest <- stRef
    nameTest <- deparse(substitute(stRef))
  } else {
    nameTest <- deparse(substitute(stTest))
  }
  nameRef <- deparse(substitute(stRef))
  

  ## stRef and stTest should be spikeTrain or repeatedTrain objects
  if (!is.spikeTrain(stRef) && !is.repeatedTrain(stRef))
    stop(paste(nameRef,
               "should be a spikeTrain or a repeatedTrain object.")
         )
  if (!is.spikeTrain(stTest) && !is.repeatedTrain(stTest))
    stop(paste(nameTest,
               "should be a spikeTrain or a repeatedTrain object.")
         )

  if (is.spikeTrain(stRef)) {
    nbRepeat <- 1
    stRef <- as.repeatedTrain(list(stRef))
    if (!is.spikeTrain(stTest)) {
      stop(paste(nameRef,"and",nameTest,
               "should be spikeTrain objects.")
         )
    } else {
      stTest <- as.repeatedTrain(list(stTest))
    }
  } else {
    nbRepeat <- length(stRef)
    if (nbRepeat != length(stTest))
      stop(paste(nameRef,"and",nameTest,
                 "should be repeatedTrain objects with identical number of repeats.")
           )
  } 
  
  isi <- function(l) unlist(lapply(l,diff))
  if (missing(laglim)) laglim <- c(-3,3)*sd(isi(stRef))

  if (missing(acquisitionWindow)) {
    acquisitionWindow <- numeric(2)
    acquisitionWindow[1] <- floor(min(c(unlist(stRef),unlist(stTest))))
    acquisitionWindow[2] <- ceiling(max(c(unlist(stRef),unlist(stTest))))
  }
  
  lagList <- function(ref,test,laglim) {
    totalRef <- 0
    totalTest <- 0
    result <- lapply(1:nbRepeat,
                     function(repIdx)
                     {
                       spikeTrainRef <- ref[[repIdx]]
                       goodRef <- ((spikeTrainRef-acquisitionWindow[1]) >= -laglim[1]) &
                       ((acquisitionWindow[2]-spikeTrainRef) >= laglim[2])
                       spikeTrainRef <- spikeTrainRef[goodRef]
                       spikeTrainTest <- test[[repIdx]]
                       ## goodTest <- ((spikeTrainTest-acquisitionWindow[1]) >= -laglim[1]) &
                       ## ((acquisitionWindow[2]-spikeTrainTest) >= laglim[2])
                       nbRefSpikes <- length(spikeTrainRef)
                       totalRef <<- totalRef + nbRefSpikes
                       ## totalTest <<- totalTest + sum(goodTest)
                       totalTest <<- totalTest + length(spikeTrainTest)
                       repResult <- lapply(1:nbRefSpikes,
                                           function(idx) {
                                             newZero <- spikeTrainRef[idx]
                                             shiftedTest <- spikeTrainTest-newZero
                                             goodTest <- laglim[1] <= shiftedTest &
                                             shiftedTest <= laglim[2]
                                             shiftedTest <- shiftedTest[goodTest]
                                             return(shiftedTest)
                                           }
                                           )
                       list(refTime=spikeTrainRef,
                            repIdx=repIdx,
                            crossTime=repResult)
                     }
                     )
    output <- vector("list",totalRef)
    gblIdx <- 1
    repIdx <- 1
    for (idx in 1:nbRepeat) {
      nbRefSpikes <- length(result[[idx]]$refTime)
      for (spikeIdx in 1:nbRefSpikes) {
        output[[gblIdx]] <- list(refTime=result[[idx]]$refTime[spikeIdx],
                                 repIdx=result[[idx]]$repIdx,
                                 crossTime=result[[idx]]$crossTime[[spikeIdx]]
                                 )
        gblIdx <- gblIdx+1
      } ## End of for loop on spikeIdx
    } ## End of for loop on idx
    output <- c(output,list(nbTestSpikes=totalTest))
    return(output)
  }

  shifted <- lagList(stRef,stTest,laglim)
  result <- list(shiftedT=shifted[-length(shifted)],
                 nbRefSpikes=length(shifted)-1,
                 nbTestSpikes=shifted[[length(shifted)]],
                 laglim=laglim,
                 acquisitionWindow=acquisitionWindow,
                 obsTime=(diff(acquisitionWindow)-diff(laglim))*nbRepeat,
                 call=match.call()
                 )
  class(result) <- "lockedTrain"
  return(result)

}

print.lockedTrain <- function(x,...) plot(x)

plot.lockedTrain <- function(x,
                             keepTime=FALSE,
                             stimTimeCourse=NULL,
                             colStim="grey80",
                             xlim,
                             pch,
                             xlab,
                             ylab,
                             main,
                             ...
                             ) {

  ## if stimTimeCourse is not NULL it should be a vector with 2 elements
  if (!is.null(stimTimeCourse)) {
    if (length(stimTimeCourse) != 2)
      stop(paste(deparse(substitute(stimTimeCourse)),
                 "should be a vector with 2 elements.")
           )
  } ## End of conditional on !is.null(stimTimeCourse)

  if (missing(xlab)) xlab <- "Lag (s)"

  nameRef <- deparse(x$call[["stRef"]])
  if (missing(ylab)) {
    if (keepTime) {
      ylab <- paste("Time of",
                    nameRef,
                    "spikes (s)")
    } else {
      ylab <- paste("Index of",
                    nameRef,
                    "spikes")
    }
  }

  if (missing(main)) {
    if (!is.null(x$call[["stTest"]])) {
      nameTest <- deparse(x$call[["stTest"]])
      main <- paste("Cross raster plot of",
                    nameRef,
                    "and",
                    nameTest
                    )
    } else {
      main <- paste("Cross raster plot of",
                    nameRef,
                    "with itself."
                    )
    } ## End of conditional on !is.null(x$call[["stTest"]])
  } ## End of conditional on missing(main)

  nbRefSpikes <- x$nbRefSpikes
  allT <- unlist(lapply(x$shiftedT, function(l) l$refTime))
  yIdx <- sort.int(allT,index.return=TRUE)$ix
  lastRefSpike <- max(allT)
  laglim <- x$laglim
  if (missing(xlim)) xlim <- laglim
  if (missing(pch)) pch <- "*"
  
  plot(laglim,
       if (keepTime) {c(0,lastRefSpike)} else {c(1,nbRefSpikes+1)},
       type="n",
       xlab=xlab,
       ylab=ylab,
       xlim=xlim,
       ylim=if (keepTime) {c(0,lastRefSpike)} else {c(1,nbRefSpikes+1)},
       bty="n",
       main=main
       )

  sortedRefSpikes <- allT[yIdx]
  if (keepTime) {
    if (!is.null(stimTimeCourse))
      rect(laglim[1],stimTimeCourse[1],
           laglim[2],stimTimeCourse[2],
           col="gray80",
           lty=0)
         
    Y <- sortedRefSpikes
  } else {
    if (!is.null(stimTimeCourse)) {
      goodOnes <- stimTimeCourse[1] <= sortedRefSpikes &
      sortedRefSpikes <= stimTimeCourse[2]
      rect(laglim[1],min((1:nbRefSpikes)[goodOnes]),
           laglim[2],max((1:nbRefSpikes)[goodOnes]),
           col="gray80",
           lty=0)
    }
    Y <- 1:nbRefSpikes
  }

  abline(v=0,col="gray30",lty=2)
  
  invisible(sapply(yIdx,
                   function(idx) {
                     xx <- x$shiftedT[[idx]]$crossTime
                     yy <- numeric(length(xx)) + Y[idx]
                     points(xx,yy,pch=pch)
                   }
                   )
            )
  
}


hist.lockedTrain <- function(x,
                             bw,
                             breaks=NULL,
                             plot=TRUE,
                             ...
                             ) {
  
  if (missing(bw)) {
    testISI <- unlist(lapply(x$shiftedT,
                             function(l) {
                               testT <- l$crossTime
                               c(diff(testT[testT<0]),
                                 diff(testT[testT>0]))
                             }
                             )
                      )
    bw <- sd(testISI)/5
  }
  
  laglim <- x$laglim
  nbRefSpikes <- x$nbRefSpikes
  allT <- unlist(lapply(x$shiftedT, function(l) l$crossTime))
  if (is.null(breaks)) {
    if (laglim[1] < 0) {
      ##l <- ((laglim[1]+bw/2)%/%bw)*bw-bw/2
      l <- laglim[1]
      ##r <- ((laglim[2]-bw/2)%/%bw)*bw+bw/2
      r <- laglim[2]
    } else {
      l <- 0
      ##r <- (laglim[2] %/% bw)*bw
      r <- laglim[2]
    }
    breaks <- seq(l,r,by=bw)
  } 
  l <- min(breaks)
  r <- max(breaks)
 
  
  allT <- allT[l < allT & allT < r]
  myHist <- hist(allT,breaks=breaks,plot=FALSE)

  nRef <- nbRefSpikes
  refFreq <- nbRefSpikes/x$obsTime
  testFreq <- x$nbTestSpikes/x$obsTime
  
  density <- myHist$counts/nRef/diff(breaks)
  if (!is.null(x$call[["stTest"]])) CCH <- TRUE
  else CCH <- FALSE

  result <- list(density=density,
                 breaks=myHist$breaks,
                 mids=myHist$mids,
                 bw=bw,
                 nRef=nRef,
                 refFreq=refFreq,
                 testFreq=testFreq,
                 obsTime=x$obsTime,
                 CCH=CCH,
                 call=match.call()
                 )
  class(result) <- "hist.lockedTrain"
  
  if (!plot) return(result)
  else plot(result,...)
  
}

plot.hist.lockedTrain <- function(x,
                                  style=c("Ogata","Brillinger"),
                                  CI,
                                  unit="s",
                                  xlab,
                                  ylab,
                                  xlim,
                                  ylim,
                                  type,
                                  pch,
                                  ...) {

  ## check CI
  ## make sure that CI is at most of length 2 otherwise
  ## keep the first 2 components
  if (missing(CI)) CI <- 0.95
  if (length(CI) > 2) CI <- CI[1:2]
  ## Check that each component of CI is between 0 and 1
  if (any(CI>=1 | CI<=0))
    stop(paste(deparse(substitute(CI)),
               "components should be in (0,1)")
         )
  bw <- x$bw
  refFreq <- x$refFreq
  testFreq <- x$testFreq
  obsTime <- x$obsTime
  lambda <- refFreq*testFreq*bw*obsTime
  lwr <- sapply(CI, function(ci) qpois((1-ci)/2,lambda))
  upr <- sapply(CI, function(ci) qpois(1-(1-ci)/2,lambda))
  lwr <- lwr/bw/obsTime/refFreq
  upr <- upr/bw/obsTime/refFreq
  if (missing(xlab)) xlab <- paste("Time Lag (",unit,")",sep="")
  if (missing(xlim)) xlim <- range(x$breaks)
  if (missing(type)) type <- "o"
  if (missing(pch)) pch <- 8
  xx <- x$mids
  CCH <- x$CCH

  if (style[1] == "Ogata") {
    yy <- x$density
    if (missing(ylab) && CCH) ylab <- "CCH"
    if (missing(ylab) && !CCH) ylab <- "ACH"
    if (missing(ylim)) {
      ylim <- c(0.9*min(c(yy,lwr)),max(c(yy,upr))*1.1)
    }

    plot(xx,
         yy,
         type="n",
         xlab=xlab,
         ylab=ylab,
         xlim=xlim,
         ylim=ylim,
         xaxs="i",
         yaxs="i")
    sapply(lwr, function(l) abline(h=l,lty=2))
    sapply(upr, function(u) abline(h=u,lty=2))       
    abline(h=testFreq)
    abline(v=0)
    lines(xx,
          yy,
          type=type,
          pch=pch,
          ...)
  } else {
    yy <- sqrt(x$density)
    if (missing(ylab) && CCH) ylab <- expression(sqrt(CCH))
    if (missing(ylab) && !CCH) ylab <- expression(sqrt(ACH))
    lwr <- sapply(lwr, function(ci) sqrt(ci))
    upr <- sapply(upr, function(ci) sqrt(ci))
    
    if (missing(ylim)) {
      ylim <- c(0,max(c(max(upr),max(yy)))*1.01)
    }

    plot(xx,
         yy,
         type="n",
         xlab=xlab,
         ylab=ylab,
         xlim=xlim,
         ylim=ylim,
         xaxs="i",
         yaxs="i")
    sapply(lwr, function(l) abline(h=l,lty=2))
    sapply(upr, function(u) abline(h=u,lty=2))       
    abline(h=sqrt(testFreq))
    lines(xx,
          yy,
          type=type,
          pch=pch,
          ...)
    
  } ## End of conditional on style
  
}

