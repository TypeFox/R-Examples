 growthTraject <- function(K, Linf, dat, lentag, lenrec, timelib, subsets=NULL,  
                          main = "Growth trajectories & fitted curve",
                          cex.lab=1.5, cex.axis=1.5, cex.main=1,
                          xlab="Relative age, yr", ylab="Length, cm",
                          xlim=NULL, ylim=NULL,ltytraject=1, lwdtraject=1,
                          coltraject=1, ltyvonB=1, lwdvonB=2, colvonB="red", 
                          returnvec=FALSE, returnlimits=FALSE, warn=TRUE, ...) {
  
  stopifnot(is.numeric(K) && K > 0, is.numeric(Linf) && Linf > 0)
  if(!missing(dat)) {
    stopifnot(dim(dat)[2] >= 3)
    if(all(c("lentag", "lenrec", "timelib") %in% names(dat))) {
      lentag <- dat$lentag
      lenrec <- dat$lenrec
      timelib <- dat$timelib
    }
    else {
      lentag <- dat[, 1]
      lenrec <- dat[, 2]
      timelib <- dat[, 3]
      if(warn == TRUE) cat("assuming lentag, lenrec & timelib in cols 1, 2, & 3, respectively\n")
    }
  }
  else {
    stopifnot(!missing(lentag) && !missing(lenrec)
              && !missing(timelib))
  }
  
  stopifnot(is.numeric(lentag) && all(lentag > 0), 
            is.numeric(lenrec) && all(lenrec > 0),
            is.numeric(timelib) && all(timelib > 0))
  if(mean(lenrec) < mean(lentag)) stop("mean recaptured length < mean length at tagging")
  
  agetag <- -log(1-lentag/Linf)/K
  agerec <- agetag + timelib
  if(missing(xlim)) xlim <- c(0, max(agerec))
  if(missing(ylim)) ylim <- c(0, max(Linf, lenrec))
  
  par(las=1, mar=c(6,6,4,2)+0.1,mgp=c(4,1,0))
  plot(2,2,xlim = xlim, ylim = ylim, xlab=xlab, 
       ylab=ylab, col=0, cex.lab=cex.lab, cex.axis=cex.axis, main=main,
       ...)
  if(!is.null(subsets)) {
    if(is.factor(subsets)) subsets <- as.integer(subsets)
    if(!is.numeric(subsets)) stop("subsets is not a factor or numeric")
    numgrp <- length(unique(subsets))
    if(length(coltraject) < numgrp && length(ltytraject) < numgrp 
       && length(lwdtraject) < numgrp) stop("Separate lines specified for each subset
                                            but not enough styles defined")
    if(length(coltraject) < numgrp) coltraject <- rep(coltraject[1], numgrp)
    if(length(ltytraject) < numgrp) ltytraject <- rep(ltytraject[1], numgrp)
    if(length(lwdtraject) < numgrp) lwdtraject <- rep(lwdtraject[1], numgrp)
    coltraject <- coltraject[subsets]
    ltytraject <- ltytraject[subsets]
    lwdtraject <- lwdtraject[subsets]
  }
  arrows(x0=agetag, y0=lentag, x1=agerec, y1=lenrec, code=0,
         lty=ltytraject, lwd=lwdtraject, col=coltraject)
  vonB <- function(x) Linf*(1-exp(-K*x))
  curve(vonB, 0, max(agetag), add=TRUE, lty=ltyvonB, lwd=lwdvonB, col=colvonB)
  
  if(returnvec == TRUE || returnlimits == TRUE) {
    answer <-list(NULL)
    if(returnvec == TRUE) {
      answer[[1]] <- data.frame(agetag=agetag, lentag=lentag, agerec=agerec, 
                                lenrec=lenrec)
      if(returnlimits == TRUE) {
        answer[[2]] <- xlim
        answer[[3]] <- ylim
      }
    }
    else {
      answer[[1]] <- xlim
      answer[[2]] <- ylim      
    }
    return(answer)
  }
}
