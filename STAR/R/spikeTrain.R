as.spikeTrain <- function(x) {
  if (inherits(x,"CountingProcessSamplePath")) x <- spikeTimes(x)
  if (!is.numeric(x)) 
    x <- as.numeric(x)
  if (!identical(length(unique(x)), length(x))) 
    stop(paste("The elements of", deparse(substitute(x)), 
               "should be all different."))
  isi <- diff(x)
  if (any(isi <= 0)) 
    stop(paste(deparse(substitute(x)), "should have strictly incresing elements."))

  class(x) <- "spikeTrain"
  x
}

is.spikeTrain <- function(x) {
  if (!("spikeTrain" %in% class(x))) return(FALSE) 
  if (!is.numeric(x)) 
    x <- as.numeric(x)
  if (!identical(length(unique(x)), length(x))) return(FALSE)
  isi <- diff(x)
  if (any(isi <= 0)) return(FALSE)
  return(TRUE)
}

plot.spikeTrain <- function(x,
                            xlab="Time (s)",
                            ylab="Cumulative Number of Events",
                            main=paste("Counting Process of", deparse(substitute(x))),
                            xlim=c(floor(x[1]), ceiling(x[length(x)])),
                            ylim=c(0, length(x) + 1),
                            do.points=ifelse(length(x)<100,TRUE,FALSE),
                            addMeanRate=TRUE,
                            addRug=TRUE,
                            ...) {
  
  if (!is.spikeTrain(x)) x <- as.spikeTrain(x)
  
  cp <- stepfun(x, 0:length(x))

  if (!addMeanRate) {
    plot(cp, xlab = xlab, ylab = ylab, main = main, do.points=do.points, ...)
  }
  else {
    plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab, 
         main = main, ...)
    abline(a = -xlim[1]*length(x)/diff(xlim), b = length(x)/diff(xlim), 
           col = 2)
    plot(cp,add=TRUE,do.points=do.points,...)
  }
  if (addRug) rug(x)
}

print.spikeTrain <- function(x,...) plot.spikeTrain(x,main="Spike Train Counting Process")

summary.spikeTrain <- function(object,timeUnit="s",digits=3,...) {
  spikeTrain <- object
  rm(object)
  if (!is.spikeTrain(spikeTrain)) stop("Not a proper spike train.")
  stRange <- range(spikeTrain)
  stNb <- length(spikeTrain)
  isi <- diff(spikeTrain)
  stStat1 <- c(mean(isi),sd(isi))
  stStat2 <- c(mean(log(isi)),sd(log(isi)))
  cat(paste("A spike train with ",
            stNb,
            " events, starting at: ",
            round(stRange[1],digits=digits),
            " and ending at: ",
            round(stRange[2],digits=digits),
            " (",timeUnit,").\n",
            "The mean ISI is: ",
            round(stStat1[1],digits=digits),
            " and its SD is: ",
            round(stStat1[2],digits=digits),
            " (",timeUnit,").\n",
            "The mean log(ISI) is: ",
            round(stStat2[1],digits=digits),
            " and its SD is: ",
            round(stStat2[2],digits=digits),
            "\n",
            "The shortest interval is: ",
            round(min(isi),digits=digits),
            "\n and the longest is: ",
            round(max(isi),digits=digits),
            " (",timeUnit,").\n",
            sep=""
            )
      )
            
}

diff.spikeTrain <- function(x,...) {
  class(x) <- NULL
  diff(x,...)
}

"[.spikeTrain" <- function(x,i) {
  xClass <- class(x)
  x <- unclass(x)
  newI <- seq(x)[i]
  if (all(diff(newI)==1)) {
    x <- x[newI]
    class(x) <- xClass
    return(x)
  }
  stop("Wrong i specification.")
}
