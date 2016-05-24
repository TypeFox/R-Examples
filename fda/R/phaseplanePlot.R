phaseplanePlot <- function(evalarg, fdobj, Lfdobj1=1, Lfdobj2=2,
      lty=c("longdash", "solid"),
      labels=list(evalarg=seq(evalarg[1], max(evalarg), length=13),
             labels=fda::monthLetters),
      abline=list(h=0, v=0, lty=2),
      xlab="Velocity", ylab="Acceleration",
                       returnMatrix=FALSE, ... ){
##
## 1.  Check 'evalarg'
##
  if(missing(evalarg))
    evalarg <- fdobj$basis$rangeval
  if(length(evalarg)<3){
    if(length(evalarg)<2)evalarg[2] <- evalarg+1
    evalarg <- seq(evalarg[1], evalarg[2], length=181)
  }
##
## 2.  Compute points to plot
##
  Eval <- sort(unique(c(evalarg, labels$evalarg)))
  D1 <- eval.fd(Eval, fdobj, Lfdobj1, returnMatrix)
  D2 <- eval.fd(Eval, fdobj, Lfdobj2, returnMatrix)
#
  nT <- length(Eval)
  n2 <- ceiling(nT/2)
##
## 3.  Set up the plot
##
  plot(range(D1), range(D2), xlab=xlab, ylab=ylab,
       type="n", ...)
  if(!is.null(abline))do.call("abline", abline)
##
## 4.  Plot the lines
##
  lines(D1[1:n2], D2[1:n2], lty=lty[1])
  lines(D1[n2:nT], D2[n2:nT], lty=lty[2])
##
## 5. Label midmonths
##
  D1. <- eval.fd(labels$evalarg, fdobj, Lfdobj1, returnMatrix)
  D2. <- eval.fd(labels$evalarg, fdobj, Lfdobj2, returnMatrix)
  text(D1., D2., labels$labels)
##
## 6.  Done
##
  out <- cbind(D1, D2)
  fd.name <- deparse(substitute(names))
  D1.name <- {
    if(is.numeric(Lfdobj1) && Lfdobj1==1)
      "Velocity"
    else
      paste(fd.name, deparse(substitute(Lfdobj1)), sep=".")
  }
  D2.name <- {
    if(is.numeric(Lfdobj2) && Lfdobj1==2)
      "Acceleration"
    else
      paste(fd.name, deparse(substitute(Lfdobj1)), sep=".")
  }
  dimnames(out) <- list(names(evalarg), c(D1.name, D2.name))
#
  invisible(cbind(D1, D2))
}
