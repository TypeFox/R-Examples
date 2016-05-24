`PacfPlot` <-
function(z, lag.max=15,...){
ylab<-"pacf"
lags<-1:lag.max
zeta<-ARToPacf((ar.burg(z,aic=FALSE,order.max=lag.max))$ar)
glim=1
plot(lags, zeta,  type="n", ylim=c(-glim,glim), xlab="lag", ylab=ylab, ...)
lines(c(0, length(zeta)), c(0,0), col="magenta")
sdZeta<-sqrt(diag(solve(InformationMatrixARz(zeta,lags)))/length(z))
segments(lags,zeta-1.96*sdZeta, lags, zeta+1.96*sdZeta, lwd=3, col="blue")
title(sub="95% confidence intervals for pacf") 
}

