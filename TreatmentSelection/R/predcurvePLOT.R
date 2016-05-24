predcurvePLOT <-
function(x, ci, ci.bounds, get.F, fixed.values, conf.bands, rho, trt.names, xlab, ylab, xlim, ylim, main, offset = .01, mar,...){ 

  fittedrisk.t0 <- x$derived.data$fittedrisk.t0
  fittedrisk.t1 <- x$derived.data$fittedrisk.t1
  marker <- x$derived.data$marker
  event <- x$derived.data$event
  trt <- x$derived.data$trt
  

  n = length(fittedrisk.t1)
  old.mar <- par()$mar

  if(is.null(mar)) mar = c(5.1, 4.1, 4.1, 9)
  par(mar=mar)  #mar=c(6.5, 4.5, 4.1, 2.1), oma=c(1.5,1,1.5,1),

   if(is.null(xlab)) xlab <- "% population below marker value"
   if(is.null(ylab)) ylab <- "risk given marker"
   if(is.null(xlim)) xlim <- c(0,1)
   if(is.null(ylim)) ylim <- c(0,1)
   if(is.null(main)) main <- "Risk curves by treatment"

    plot(NULL, 
          ylab = ylab,
          xlab = xlab,
          xlim = xlim, 
          ylim = ylim,
          type = "n", 
          main = main, ...)

 legend(x=xlim[2]+diff(xlim)/15, y = quantile(ylim, prob = .75), legend = trt.names, lty = c(2, 1),lwd=c(2,2), bty="n", cex = 1, xpd = TRUE)
  
  #browser()
  F.Y <- get.F(marker, event, trt, rho = rho)

  stepF.Y <- rep(sort(F.Y), c(rep(2, n-1),1))

  x.t0 <- stepF.Y
  x.t1 <- stepF.Y  

  y.t0 <- fittedrisk.t0[rep(order(F.Y),c(1, rep(2, n-1)))]
  y.t1 <- fittedrisk.t1[rep(order(F.Y),c(1, rep(2, n-1)))]
  if(!is.null(ci.bounds)){
  ci.bounds <- matrix(ci.bounds, ncol=length(fixed.values), nrow = 4)
  if(substr(ci, 1,1)=="h"){
  #the indices of fixed values that fall between min(fittedrisk.t0) and max(fittedrisk.t0) ...same for t1
  index.fix.t0   <- (fixed.values<= max(fittedrisk.t0) & fixed.values >= min(fittedrisk.t0)) 
  index.fix.t1   <- (fixed.values<= max(fittedrisk.t1) & fixed.values >= min(fittedrisk.t1)) 
  }else{

  index.fix.t0   <- (fixed.values<= max(F.Y) & fixed.values >= min(F.Y)) 
  index.fix.t1   <- (fixed.values<= max(F.Y) & fixed.values >= min(F.Y)) 
  }
  
  shade(ci.bounds[1:2,index.fix.t0], fixed.values[index.fix.t0], type = substr(ci, 1, 1), bands = conf.bands)
  shade(ci.bounds[3:4,index.fix.t1], fixed.values[index.fix.t1] +offset, type = substr(ci, 1, 1), bands = conf.bands, lty = 2)
  }
  lines(x.t0, y.t0, type = "l", lwd=2)  
  lines(x.t1, y.t1, lty=2, lwd = 2) 

  #axis(1,line=4,at=c(0,.2,.4,.6,.8,1),
  #       round(quantile(marker,probs=c(0,.2,.4,.6,.8,1)),2))
  #mtext("Marker Percentile = F(marker)", 1,  line = 2,  cex=1.2)
  #mtext("Marker Value = marker",         1,  line = 6,  cex=1.2)
  par(mar=old.mar)
  return(cbind(x.t0, y.t0, x.t1, y.t1))
}
