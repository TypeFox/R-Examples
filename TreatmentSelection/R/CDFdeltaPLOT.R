CDFdeltaPLOT <-
function(x, ci, ci.bounds, get.F, fixed.values,conf.bands,  rho, xlab, ylab, xlim, ylim, main, mar,  ...){ 

  trt.effect <- x$derived.data$trt.effect
  marker <- x$derived.data$marker
  event <- x$derived.data$event
  trt <- x$derived.data$trt
  n = length(trt.effect)


    cen <- mean(c(min(trt.effect, na.rm=TRUE), max(trt.effect, na.rm=TRUE)))
    ran <- max(trt.effect, na.rm=TRUE) - min(trt.effect, na.rm=TRUE)

    if(substr(ci, 1,1) =="h"){ cen <- mean(c(min(c(trt.effect, ci.bounds), na.rm=TRUE), max(c(trt.effect, ci.bounds), na.rm=TRUE))); ran <- max(c(trt.effect, ci.bounds), na.rm=TRUE) - min(c(trt.effect, ci.bounds), na.rm=TRUE)}
    

    ran <- ran*1.1


   old.mar = par()$mar
   if(is.null(mar)) mar = c(5.1, 4.1, 4.1, 2.1)
   par(mar=mar)  #mar=c(6.5, 4.5, 4.1, 2.1), oma=c(1.5,1,1.5,1),

    mylim <- c(cen-ran/2, cen+ran/2)
  
   if(is.null(xlab)) xlab <- "treatment effect"
   if(is.null(ylab)) ylab <- "% population below treatment effect"
   if(is.null(xlim)) xlim <- mylim
   if(is.null(ylim)) ylim <- c(0,1)
   if(is.null(main)) main <- "Treatment effect distribution"

    plot(NULL, 
          ylab = ylab,
          xlab = xlab,
          xlim = xlim, 
          ylim = ylim,
          type = "n", 
          main = main, ...)

  F.D <- get.F(trt.effect, event, trt, rho = rho)

  
  if(!is.null(ci.bounds)){ 
  ci.bounds <- matrix(ci.bounds, ncol=length(fixed.values), nrow = 2)

  if(substr(ci, 1,1)=="h"){
    index.fix  <- (fixed.values<= max(F.D) & fixed.values >= min(F.D)) 
    width = 5
  }else{
    index.fix  <- (fixed.values<= max(trt.effect) & fixed.values >= min(trt.effect)) 
    width = .05
  }
  
  shade(ci.bounds[,index.fix], fixed.values[index.fix], type = substr(ci, 1, 1), bands = conf.bands)
  }

  
  stepF.D <- c(rep(F.D[order(trt.effect)], c(rep(2, n-1), 1)))         
 
  myy <- stepF.D  
  myx <- trt.effect[rep(order(trt.effect),c(1, rep(2, n-1)))]
  
  lines(myx, myy,type = "l", lwd=2)  

  par(mar = old.mar)

  return(cbind("x" =myx, "y" =myy))
}
