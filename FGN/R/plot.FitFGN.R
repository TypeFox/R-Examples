`plot.FitFGN` <-
function(x, maxLag = 30, ...){
res <- x$res
n <- length(res)
if (maxLag >= n)
	maxLag <- ceiling(n/4)
ra <- (acf(res, plot=FALSE, lag.max=maxLag)$acf)[-1]
MXL <- length(ra)
layout(matrix(c(1,2,1,2),ncol=2))
clim<-1.96
sdra<-sqrt(1/n)
glim <- max(abs(ra),clim*sdra)*1.05
lags <- 1:MXL
plot(lags, ra, type="h", ylim=c(-glim,glim), xlab="lag", ylab="Acf", 
    main="Residual Acf and 5% Significance Limit")
lines(c(0, length(ra)), c(0,0), col="magenta")
abline(h=clim*sdra, col="blue")
abline(h=-clim*sdra, col="blue")
#title(main=x$DataTitle)
#Ljung-Box
Q <- x$LjungBoxQ
m <- Q[,1]
p <- Q[,3]
plot(m,p,ylim=c(0,1),xlab="lag", ylab="p-value", main="Ljung-Box Portmanteau Test")
abline(h=0.05, lty=2, col=3)
}
