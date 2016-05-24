isattest <-
function(x, hnull=0, lr=FALSE, ci.pval = 0.99, plot=TRUE, plot.turn = FALSE, biascorr=FALSE){
  
  trend.incl <- FALSE
  if(!is.null(as.list(x$call)$tis)){
    if (as.list(x$call)$tis) {
      stop("isat.test currently not implemented for trend-indicator saturation")
      trend.incl <- TRUE
    }
  }
  
  arcall <- as.list(x$call)$ar
  x.var <- isatvar(x,lr=lr)
  
  if (biascorr==TRUE){
    
    if (!is.null(as.list(x$call)$mxreg) | !is.null(arcall) | trend.incl){
      
      biascorr <- FALSE
      print("Bias Correction not applicable in isat regression with additional non-step covariates. Has been set to FALSE.")
    }
  }
  
  
  T <- dim(x$aux$mX)[1]
  N <- dim(x$aux$mX)[2]
  
  crit <- abs(qt((1-ci.pval)/2, T-N))
  bias.low <- matrix(0, T, 1)
  bias.high <- matrix(0, T, 1)
  
  ci.low <- matrix(0, T, 1)
  ci.high <- matrix(0, T, 1)
  x.mean <- matrix(0, T, 1)
  
  if (lr == TRUE & !is.null(arcall))
  {
    x.is.lr <- x.var$lr.path
    x.is.const <- x.var$const.path
    
  } else {
    
    x.is.lr <- NA
    
    if (biascorr){
      
      #is2$aux$t.pval
      
      xbias <-biascorr(b=x.var$const.path, b.se=x.var$const.se, p.alpha = x$aux$t.pval, T=length(x.var$const.path))
      
      #xbias <-biascorr(b=x.var$const.path, b.se=x.var$const.se, p.alpha = as.list(x$call)$t.pval, T=length(x.var$const.path))
      x.is.const <- xbias$beta.2step
      
    } else {
      x.is.const <- x.var$const.path
      
    }
    
    
  }
  
  
  if (lr == TRUE & !is.null(arcall))
  {
    
    ci.low <- x.var$lr.path-crit*x.var$lr.se
    ci.high <- x.var$lr.path+crit*x.var$lr.se
    
    bias.low[which((ci.low) > hnull)] <- 1
    bias.high[which((ci.high) < hnull)] <- 1
    
    bias.low <- bias.low*(x.var$lr.path-hnull)
    bias.high <- bias.high*(x.var$lr.path-hnull)
    
    x.mean <- x.var$lr.path
    
  } else {
    
    ci.low <- x.is.const-crit*x.var$const.se
    ci.high <- x.is.const+crit*x.var$const.se
    
    bias.low[which((ci.low) > hnull)] <- 1
    bias.high[which((ci.high) < hnull)] <- 1
    
    
    bias.low <- bias.low*(x.is.const-hnull)
    bias.high <- bias.high*(x.is.const-hnull)
    
    x.mean <- x.is.const
    
  }
  
  ###determining the turning points
  time <- x$aux$y.index
  bias.sum.ar <- bias.low+bias.high
  lr.path.d <- diff(bias.sum.ar)
  
  
  if(all(lr.path.d==0)){
    plot.turn <- FALSE
    turn.ar <- NULL
  } else {
    turn.ar <- time[which(lr.path.d != 0)]+1
  }
  
  turn.ar.y <- bias.sum.ar[which(lr.path.d != 0)]
  
  turn.x.lab <- turn.ar
  turn.x <- turn.ar
  
  fitted <- x$mean.fit
  actual <- zoo(x$aux$y, order.by=x$aux$y.index)
  
  ylabel_a <- "Series"
  ylabel_b <- "Bias"
  
  
  par(mfrow=c(2,1), mar = c(2, 4,1,3))
  Ylim_main <- c(min(actual, na.rm=TRUE)*1.2,max(actual, na.rm=TRUE)*1.2)
  Ylim_bias <- c(min(bias.high, na.rm=TRUE)*1.2,max(bias.low, na.rm=TRUE)*1.2)
  
  if (plot){
    
    plot(time, x.mean, ylim=Ylim_main, col="blue", title(main=NULL, xlab=NULL), xlab=NA, ylab=ylabel_a, sub=NA, type="l")
    lines(ci.low, col="blue", lty=2)
    lines(ci.high, col="blue", lty=2)
    lines(actual)
    abline(a =hnull, b=0, col="black", lty=3, lwd=2)
    
    plot(time, bias.low, type="h", col="red", ylim=Ylim_bias, title(main=NULL, xlab=NULL), xlab=NA, ylab=ylabel_b, sub=NA)
    lines(bias.high, type="h", col="blue")
    
    if ( plot.turn ){
      text(turn.x.lab, y=turn.ar.y, x=turn.x, pos=4, offset=-0.5, cex=0.8)
    }
    
  }
  
  
  if (lr==TRUE & !is.null(arcall))
  {
    mean.var <- cbind(ci.low, ci.high, bias.low,  bias.high)
    colnames(mean.var) <- c("ci.low", "ci.high", "bias.high",  "bias.low")
    
    
  } else {
    mean.var <- cbind( ci.low, ci.high, bias.low,  bias.high)
    colnames(mean.var) <- c("ci.low", "ci.high", "bias.high",  "bias.low")
    
  }
  
  
  return(mean.var)
  
}
