HWQqplot <- function(X,nsim=100,fit="curve",logplot=FALSE,main="Q-Q plot for HWE",mm=NULL,pvaluetype="selome",...) {
pvals <- HWExactMat(X,pvaluetype=pvaluetype)$pvalvec    
out <- sort(pvals,index.return=TRUE)
pvals <- out$x
X <- X[out$ix,]
n <- length(pvals)
epvals <- getpvals(X,nsim)
a <- quantile(pvals,probs=c(0.25,0.75))
# there are no theoretical reference quantiles
if (!logplot) {
  plot(pvals,pvals,pch=19,xlab="Expected p value",ylab="Observed p value",xlim = c(0,1), ylim = c(0,1) , type="n", main = main, ...)
  for(i in 1:nsim) {
    b <- quantile(epvals[,i],probs=c(0.25,0.75))
    slope <- diff(a)/diff(b)
    int <- a[1] - slope * b[1]
    #  mcat(int,slope)
    if(fit=="line") abline(int, slope, lwd = 1, col = "grey")
    if(fit=="curve") points(sort(epvals[,i]),pvals,type="l",col="grey")
  }
  abline(0,1,col="green",lwd=2)
}
  else {
    # plot in logarithmic (base 10) scale.
   lpvals <- -log10(pvals)
   lepvals <- -log10(epvals)        
   if(is.null(mm)) mm <- max(c(lpvals, lepvals))
   plot(lpvals, lpvals, xlab = "-log(Expected p value)", 
   ylab = "-log(Observed p value)", xlim = c(0, mm), 
   ylim = c(0, mm), main = main, type="n", ...)
   for(i in 1:nsim) {
       qp <- quantile(lpvals, c(0.25, 0.75))
       qe <- quantile(lepvals[,i], c(0.25, 0.75))
       slope <- diff(qp)/diff(qe)
       int <- qp[1] - slope * qe[1]
       if(fit=="line") abline(int, slope, lwd = 1, col = "grey")
       if(fit=="curve") points(sort(lepvals[,i]),sort(lpvals),type="l",col="grey") 
   }
    abline(0,1,col="green",lwd=2)
  }
}
