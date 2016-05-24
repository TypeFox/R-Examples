trteffectPLOT <-
function(x, ci, ci.bounds, get.F, fixed.values, conf.bands,  rho, xlab, ylab, xlim, ylim, main, markerTWO=FALSE, lty = 1, mar,  ...){ 
  
  trt.effect <- x$derived.data$trt.effect
  marker <- x$derived.data$marker
  event <- x$derived.data$event
  trt <- x$derived.data$trt
  n = length(trt.effect)
  


  

 if(!markerTWO){
 old.mar <- par()$mar

  if(is.null(mar)) mar = c(5.1, 4.1, 4.1, 9)
  par(mar=mar)  #mar=c(6.5, 4.5, 4.1, 2.1), oma=c(1.5,1,1.5,1),

    #browser()
    cen <- mean(c(min(trt.effect, na.rm=TRUE), max(trt.effect, na.rm=TRUE)))
    ran <- max(trt.effect, na.rm=TRUE) - min(trt.effect, na.rm=TRUE)

    if(substr(ci, 1,1) =="v"){ cen <- mean(c(min(c(trt.effect, ci.bounds), na.rm=TRUE), max(c(trt.effect, ci.bounds), na.rm=TRUE))); ran <- max(c(trt.effect, ci.bounds), na.rm=TRUE) - min(c(trt.effect, ci.bounds), na.rm=TRUE)}
    

    ran <- ran*1.1
    mylim <- c(cen-ran/2, cen+ran/2)

   if(is.null(xlab)) xlab <- "% population below treatment effect"
   if(is.null(ylab)) ylab <- "treatment effect"
   if(is.null(xlim)) xlim <- c(0,100)
   if(is.null(ylim)) ylim <- mylim
   if(is.null(main)) main <- "Treatment effect distribution"

    plot(NULL, 
          ylab = ylab,
          xlab = xlab,
          xlim = xlim, 
          ylim = ylim,
          type = "n", 
          main = main, ...)
     legend(x=xlim[2]+diff(xlim)/15, y = quantile(ylim, prob = .75), legend = c("Average", "Zero"), lty = c(3,4), col = c("black", "black"), bty="n", cex = 1, xpd = TRUE, title= "Treatment Effect")
     lines(c(0, 100), c(0,0), lty = 4, col="black")
     lines(c(0, 100), rep(mean(event[trt==0])-mean(event[trt==1]), 2), lty = 3)
    
  par(mar = old.mar)
  }

  
  F.D <- get.F(trt.effect, event, trt, rho = rho)*100

 
  if(!is.null(ci.bounds)){
 
  ci.bounds <- matrix(ci.bounds, ncol=length(fixed.values), nrow = 2)

  if(substr(ci, 1,1)=="v"){
    index.fix  <- (fixed.values<= max(F.D) & fixed.values >= min(F.D)) 
  }else{
    index.fix  <- (fixed.values<= max(trt.effect) & fixed.values >= min(trt.effect)) 
  }
  
  
  shade(ci.bounds[,index.fix], fixed.values[index.fix], type = substr(ci, 1, 1), bands = conf.bands, lty)
  }
  
  stepF.D <- c(rep(F.D[order(F.D)], c(rep(2, n-1), 1)))         
 
  myx <- stepF.D  
  myy <- trt.effect[rep(order(F.D),c(1,rep(2, n-1)))]
  lines(myx, myy,type = "l", lwd=2, lty = lty)  
  
  return(cbind("x" =myx, "y" =myy))
}
