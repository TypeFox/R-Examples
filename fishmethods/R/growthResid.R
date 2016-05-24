growthResid <- function(K, Linf, dat, lentag, lenrec, timelib, graph=1, 
                        main = "Residuals of growth increments",
                        cex.lab=1.5, cex.axis=1.5, cex.main=1,
                        xlab1="Relative age, yr", xlab2="Time at liberty, yr",
                        ylab="Observed - expected increment",
                        xlim1=NULL, xlim2=NULL, ylim=NULL, col=1, returnvec=FALSE, 
                        returnlimits=FALSE, warn=TRUE, ...) {
  xlimout<-NULL
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
  
  obs <- lenrec - lentag
  agetag <- -log(1-lentag/Linf)/K
  agerec <- agetag + timelib      
  pred <- function(lentag, agerec, K, Linf) {
    Lr <- Linf*(1 - exp(-K*agerec)) # predicted length at recapture
    predict_incr <- Lr - lentag  # predicted length increment
    return(predict_incr)
  }
  
  preds <- pred(lentag, agerec, K, Linf)
  resids <- obs - preds
  
  if(is.null(xlim1)) xlim1 <- c(0, max(agetag))
  if(is.null(xlim2)) xlim2 <- c(0, max(timelib))
  if(is.null(ylim)) ylim <- c(min(resids), max(resids))
  
  par(las=1)
  if(graph==1){
     plot(agetag,resids, pch=16, col=rgb(0,0,0,.1), xlab=xlab1, ylab=ylab, 
       main=main, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
       xlim=xlim1, ylim=ylim)
     xlimout<-xlim1
    abline(h=0, lty=3)
   }
  if(graph==2){
    plot(timelib,resids, pch=16, col=rgb(0,0,0,.1), xlab=xlab2, ylab=ylab, 
       main=main, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
       xlim=xlim2, ylim=ylim)
       xlimout<-xlim2
    abline(h=0, lty=3)
  }
  if(returnvec == TRUE || returnlimits == TRUE) {
    answer <-list(NULL)
    if(returnvec == TRUE) {
      answer[[1]] <- data.frame(agetag=agetag, resids=resids)
      names(answer)[[1]]<-c("Relative Age and Residuals")
      if(returnlimits == TRUE) {
        answer[[2]] <- xlimout
        answer[[3]] <- ylim
        if(graph==1){
          names(answer)[[2]]<-c("Graph 1 X Limits")
          names(answer)[[3]]<-c("Graph 1 Y Limits")
        } 
        if(graph==2){ 
          names(answer)[[2]]<-c("Graph 2 X Limits")
          names(answer)[[3]]<-c("Graph 2 Y Limits")
        }
      }
    }
    else {
      answer[[1]] <- xlimout
      answer[[2]] <- ylim 
      if(graph==1){
        names(answer)[[1]]<-c("Graph 1 X Limits")
        names(answer)[[2]]<-c("Graph 1 Y Limits")
      } 
      if(graph==2){ 
        names(answer)[[1]]<-c("Graph 2 X Limits")
        names(answer)[[2]]<-c("Graph 2 Y Limits")
      }
    }
    return(answer)
  }
}
