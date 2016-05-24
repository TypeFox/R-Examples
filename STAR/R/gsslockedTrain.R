.mkDiscreteLockedTrain <- function(lockedTrain,
                                   bw
                                   ) {
  ## check that lockedTrain is a "lockedTrain" object
  if (!inherits(lockedTrain,"lockedTrain"))
    stop("lockedTrain should be a \"lockedTrain\" object.")
  ## make a histogram without generating a plot
  lTh <- hist(lockedTrain,bw=bw,plot=FALSE)
  Time <- lTh$mids
  nRef <- lTh$nRef
  testFreq <- lTh$testFreq
  bwV <- diff(lTh$breaks)
  Count <- lTh$density*nRef*bwV
  result <- data.frame(Count=Count,
                       Time=Time)
  attr(result,"nRef") <- nRef
  attr(result,"testFreq") <- testFreq
  attr(result,"bwV") <- bwV
  attr(result,"CCH") <- lTh$CCH
  return(result)
}


gsslockedTrain <- function(lockedTrain,
                           bw=0.001,
                           ...)
### smooths a lockedTrain object
### uses the gss model with the Poisson family after
### binning the object.
{
  df <- .mkDiscreteLockedTrain(lockedTrain,bw)
  PoissonFit <- gssanova(Count ~ Time, family="poisson", data=df, ...)

  result <- list(gssFit=PoissonFit,
                 Time=df$Time,
                 nRef=attr(df,"nRef"),
                 testFreq=attr(df,"testFreq"),
                 bwV=attr(df,"bwV"),
                 CCH=attr(df,"CCH"),
                 call=match.call())
  
  class(result) <- "gsslockedTrain"
  result
}

gsslockedTrain0 <- function(lockedTrain,
                            bw=0.001,
                            ...)
### smooths a lockedTrain object
### uses the gss model with the Poisson family after
### binning the object.
{
  df <- .mkDiscreteLockedTrain(lockedTrain,bw)
  PoissonFit <- gssanova0(Count ~ Time, family="poisson", data=df, ...)

  result <- list(gssFit=PoissonFit,
                 Time=df$Time,
                 nRef=attr(df,"nRef"),
                 testFreq=attr(df,"testFreq"),
                 bwV=attr(df,"bwV"),
                 CCH=attr(df,"CCH"),
                 call=match.call())
  
  class(result) <- "gsslockedTrain0"
  result
}

gssObj.gsslockedTrain <- function(object, ...) {
  object$gssFit
}

gssObj.gsslockedTrain0 <- gssObj.gsslockedTrain

summary.gsslockedTrain <- function(object, ...) {
  summary(gssObj(object))
}

summary.gsslockedTrain0 <- summary.gsslockedTrain

print.gsslockedTrain <- function(x, ...) {
  print(gssObj(x))
}

print.gsslockedTrain0 <- print.gsslockedTrain

plot.gsslockedTrain <- function(x,
                                xlab, ylab, main, xlim, ylim,
                                col,lwd,
                                ...) {

  if (missing(xlab)) 
    xlab <- "Time (s)"
  if (missing(ylab)) 
    ylab <- "PDF of a Test Neuron Spike (1/s)"
  if (missing(main)) {
    nameList <- deparse(x$call[["lockedTrain"]])
    if (x$CCH) main <- paste(nameList, "smoothed cross-intensity")
    else main <- paste(nameList, "smoothed auto-intensity")
  }

  gssFitP <- predict(gssObj(x),newdata=data.frame(Time=x$Time),se.fit=TRUE)
  Y <- exp(gssFitP$fit)/x$nRef/x$bwV
  Yp <- exp((gssFitP$fit+1.96*gssFitP$se.fit))/x$nRef/x$bwV
  Ym <- exp((gssFitP$fit-1.96*gssFitP$se.fit))/x$nRef/x$bwV
  X <- x$Time

  if (missing(xlim)) xlim <- range(X)
  if (missing(ylim)) ylim <- c(min(c(Ym,x$testFreq)),max(c(Yp,x$testFreq)))
  if (missing(col)) col <- 2
  if (missing(lwd)) lwd <- 2
  
  plot(X,Y,type="n",
       xlab=xlab,ylab=ylab,
       xlim=xlim,ylim=ylim,
       main=main,
       ...)
  abline(h=x$testFreq)
  abline(v=0)
  lines(X,Yp,lty=2)
  lines(X,Ym,lty=2)
  lines(X,Y,col=col,lwd=lwd)
  
}

plot.gsslockedTrain0 <- plot.gsslockedTrain
