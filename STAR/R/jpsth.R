jsd <- function(xRT,
                yRT,
                acquisitionWindow,
                xlab,
                ylab,
                main,
                pch=".",
                ...
                )
### computes a joint peristimulus scatter diagram
### xRT and yRT should be repeatedTrain objects.
### The spike times of xRT are used to get the x coordinate
### The spike times of yRT are used to get the y coordinate
### acquisitionWindow the acquisition window used or a subset
### of it.
{
  xRTn <- deparse(substitute(xRT))
  yRTn <- deparse(substitute(yRT))
  
  ## check that xRT and yRT are repeatedTrain objects of the
  ## same length
  if (!is.repeatedTrain(xRT))
    stop(paste(xRTn,"should be a \"repeatedTrain\" object."))
  if (!is.repeatedTrain(yRT))
    stop(paste(yRTn,"should be a \"repeatedTrain\" object."))
  nbTrials <- length(xRT)
  if (nbTrials != length(yRT))
    stop(paste(xRTn,"and", yRTn,
               "should be \"repeatedTrain\" objects of the same length.")
         )

  if (missing(acquisitionWindow)) {
    acquisitionWindow <- range(c(unlist(xRT),unlist(yRT)))
    acquisitionWindow <- c(floor(acquisitionWindow[1]),
                           ceiling(acquisitionWindow[2])
                           )
  }

  if (missing(xlab))
    xlab <- paste(xRTn, "spike times (s)")

  if (missing(ylab))
    ylab <- paste(yRTn, "spike times (s)")

  if (missing(main))
    main <- paste("Joint peri-stimulus scatter diagram of:",
                  xRTn,"and",yRTn)

  plot(acquisitionWindow,
       acquisitionWindow,
       type="n",
       xlab=xlab,
       ylab=ylab,
       main=main,
       ...)

  invisible(sapply(1:nbTrials,
                   function(tIdx) {
                     y <- as.numeric(yRT[[tIdx]])
                     X <- as.numeric(xRT[[tIdx]])
                     sapply(X,
                            function(x) {
                              x <- rep(x,length(y))
                              points(x,y,pch=pch)
                            }
                            )
                   }
                   )
            )
  
}

jpsth <- function(xRT,
                  yRT,
                  xBreaks,
                  yBreaks,
                  acquisitionWindow,
                  nbEvtPerBin=50
                  )
### computes a joint peristimulus time-histogram
### xRT and yRT should be repeatedTrain objects.
### The spike times of xRT are used to get the x coordinate
### The spike times of yRT are used to get the y coordinate
### xBreaks can be a single number specifying the number
### of equally larged bins to use to split acquisitionWindow
### or a vector of boundaries between bins for the x axis.
### yBreaks can be a single number specifying the number
### of equally larged bins to use to split acquisitionWindow
### or a vector of boundaries between bins for the y axis.
### acquisitionWindow the acquisition window used or a subset
### of it.
{
  ## check that xRT and yRT are repeatedTrain objects of the
  ## same length
  if (!is.repeatedTrain(xRT))
    stop("xRT should be a \"repeatedTrain\" object.")
  if (!is.repeatedTrain(yRT))
    stop("yRT should be a \"repeatedTrain\" object.")
  nbTrials <- length(xRT)
  if (nbTrials != length(yRT))
    stop("xRT and yRT should be \"repeatedTrain\" objects of the same length.")

  if (missing(acquisitionWindow)) 
    acquisitionWindow <- range(c(unlist(xRT),unlist(yRT)))
  
  l <- floor(acquisitionWindow[1])
  r <- ceiling(acquisitionWindow[2])

  if (missing(xBreaks) && missing(yBreaks)) {
    xRT <- lapply(xRT, function(st) as.numeric(st[l < st & st < r]))
    yRT <- lapply(yRT, function(st) as.numeric(st[l < st & st < r]))

    xFreq <- sum(sapply(xRT, length))/(r-l)
    yFreq <- sum(sapply(yRT, length))/(r-l)

    bw <- sqrt(nbEvtPerBin/xFreq/xFreq)
    xBreaks <- seq(l,r,bw)
    yBreaks <- xBreaks
    r <- min(r,xBreaks[length(xBreaks)])
  } ## End of conditional on missing(xBreaks) && missing(yBreaks)

  if (length(xBreaks) == 1) {
    ## xBreaks is interpreted as a bin width
    bw <- xBreaks
    xBreaks <- seq(l,r,bw)
  } ## End of conditional on length(xBreaks) == 1
  
  if (length(yBreaks) == 1) {
    ## yBreaks is interpreted as a bin width
    bw <- yBreaks
    yBreaks <- seq(l,r,bw)
  } ## End of conditional on length(yBreaks) == 1

  if (missing(xBreaks) && !missing(yBreaks)) xBreaks <- yBreaks
  if (missing(yBreaks) && !missing(xBreaks)) yBreaks <- xBreaks

  l <- min(c(xBreaks,yBreaks))
  r <- max(c(xBreaks,yBreaks))
  xRT <- lapply(xRT, function(st) as.numeric(st[l < st & st < r]))
  yRT <- lapply(yRT, function(st) as.numeric(st[l < st & st < r]))

  xTotal <- sum(sapply(xRT, length))
  xFreq <- xTotal/(r-l)
  yTotal <- sum(sapply(yRT, length))
  yFreq <- yTotal/(r-l)

  xRT <- lapply(xRT, function(x) findInterval(x,xBreaks))
  yRT <- lapply(yRT, function(y) findInterval(y,yBreaks))

  counts <- matrix(as.integer(0),
                   nrow=length(xBreaks)-1,
                   ncol=length(yBreaks)-1)

  for (tIdx in 1:nbTrials) {
    x <- xRT[[tIdx]]
    y <- yRT[[tIdx]]
    for (i in seq(x)) counts[x[i],y] <- counts[x[i],y] + 1
  } ## End of for loop on tIdx

  density <- counts/(diff(xBreaks) %o% diff(yBreaks))/xTotal/yTotal

  xMids <- xBreaks[-length(xBreaks)]+diff(xBreaks)/2
  yMids <- yBreaks[-length(yBreaks)]+diff(yBreaks)/2

  result <- list(counts=counts,
                 density=density,
                 xMids=xMids,
                 yMids=yMids,
                 xBreaks=xBreaks,
                 yBreaks=yBreaks,
                 xTotal=xTotal,
                 yTotal=yTotal,
                 xFreq=xFreq,
                 yFreq=yFreq,
                 nbTrials=nbTrials,
                 acquisitionWindow=c(l,r),
                 call=match.call()
                 )

  class(result) <- "jpsth"
  return(result)
                 
}


contour.jpsth <- function(x,
                          xlab,
                          ylab,
                          main,
                          ...)
### contour method or jpsth objects
{
  xN <- deparse(substitute(x))
  ## check x
  if (!inherits(x,"jpsth"))
    stop(paste(xN,"should be a jpsth object."))

  if (missing(xlab))
    xlab <- paste(deparse(x$call[["xRT"]]),"spike times (s)")
  if (missing(ylab))
    ylab <- paste(deparse(x$call[["yRT"]]),"spike times (s)")
  if (missing(main))
    main <- paste("JPSTH of",
                  deparse(x$call[["xRT"]]),
                  "and",
                  deparse(x$call[["yRT"]])
                  )
  
  contour(x=x$xMids,
          y=x$yMids,
          z=x$density,
          xlab=xlab,
          ylab=ylab,
          main=main,
          ...)

}

image.jpsth <- function(x,
                        xlab,
                        ylab,
                        main,
                        ...)
### image method or jpsth objects
{
  xN <- deparse(substitute(x))
  ## check x
  if (!inherits(x,"jpsth"))
    stop(paste(xN,"should be a jpsth object."))

  if (missing(xlab))
    xlab <- paste(deparse(x$call[["xRT"]]),"spike times (s)")
  if (missing(ylab))
    ylab <- paste(deparse(x$call[["yRT"]]),"spike times (s)")
  if (missing(main))
    main <- paste("JPSTH of",
                  deparse(x$call[["xRT"]]),
                  "and",
                  deparse(x$call[["yRT"]])
                  )
  
  image(x=x$xMids,
        y=x$yMids,
        z=x$density,
        xlab=xlab,
        ylab=ylab,
        main=main,
        ...)

}

persp.jpsth <- function(x,
                        xlab,
                        ylab,
                        main,
                        ...)
### persp method or jpsth objects
{
  xN <- deparse(substitute(x))
  ## check x
  if (!inherits(x,"jpsth"))
    stop(paste(xN,"should be a jpsth object."))

  if (missing(xlab))
    xlab <- paste(deparse(x$call[["xRT"]]),"spike times (s)")
  if (missing(ylab))
    ylab <- paste(deparse(x$call[["yRT"]]),"spike times (s)")
  if (missing(main))
    main <- paste("JPSTH of",
                  deparse(x$call[["xRT"]]),
                  "and",
                  deparse(x$call[["yRT"]])
                  )
  
  persp(x=x$xMids,
        y=x$yMids,
        z=x$density,
        xlab=xlab,
        ylab=ylab,
        main=main,
        ...)

}


jpsth2df <- function(object)
### function formating a jpsth object
### into a data frame suitable for use
### in glm or gam for instance
{
  objectN <- deparse(substitute(object))
  ## check object
  if (!inherits(object,"jpsth"))
    stop(paste(objectN,"should be a jpsth object."))

  Count <- object$counts
  CountD <- dim(Count)
  dim(Count) <- NULL
  X <- rep(object$xMids,CountD[2])
  Y <- rep(object$yMids,each=CountD[1])
  result <- data.frame(Count=Count,
                       X=X,
                       Y=Y)
  
  attr(result,"xBreaks") <- object$xBreaks
  attr(result,"yBreaks") <- object$yBreaks
  attr(result,"xTotal") <- object$xTotal
  attr(result,"yTotal") <- object$yTotal
  attr(result,"nbTrials") <- object$nbTrials
  attr(result,"acquisitionWindow") <- object$acquisitionWindow
  attr(result,"originalCall") <- object$call

  return(result)
  
}
