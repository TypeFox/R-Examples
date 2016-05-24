
FDRCurve <- function(object, q, threshold=1, plot = TRUE, xlim , ylim , xlab , 
                     ylab , main , ...)  {
   ### object - an object of class "rvalues"
   
   if(missing(xlim)) {
       xlim <- range(object$aux$alpha.grid)
   }
   if(missing(ylim)) {
       ylim <- c(0,.8)
   }
   if(missing(xlab)) {
       xlab <- "r-value threshold"
   }
   if(missing(ylab)) {
       ylab <- "FDR"
   }
   if(missing(main)) {
       main <- "Estimated FDR by r-value threshold"
   }
   
   if(object$aux$prior=="conjugate") {
       jstar <- which.min(abs(object$aux$theta.quantiles - threshold))
   }
   else if(object$aux$prior=="nonparametric") {
   ###  For the nonparametric case:
       alpha.star <- 1 - object$aux$mixcdf(threshold) 
       jstar <- which.min(abs(object$aux$alpha.grid - alpha.star))
   }
    
   Vstar <- 1 - object$aux$V[,jstar]
   threshgrid <- seq(.001,.999,length.out=100) 
   nthresh <- length(threshgrid)
   fdr <- threshgrid
   for(j in 1:nthresh) {
        thresh <- (object$aux$unsorted$RValue <= threshgrid[j])
        fdr[j] <- sum(Vstar*thresh)/max(sum(thresh),1)
   }
   exceed_ind <- which.min(abs(fdr - q))

   ans <- list()
   ans$fdrcurve <- approxfun(threshgrid,fdr)
   ans$Rval.cut <- threshgrid[exceed_ind]
   if(plot)  {
      plot(threshgrid,fdr,type="n",xlim=xlim,ylim=ylim,xlab=xlab,
           ylab=ylab,main=main, ...)
      lines(threshgrid,fdr,lwd=2)
      abline(q,0,lwd=2,lty=2)
   
      invisible(ans)
   }
   else {
      return(ans)
   }
}
