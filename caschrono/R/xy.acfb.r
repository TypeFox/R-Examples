xy.acfb=function(y, lag.max=40, numer=TRUE)
{
num=length(y)
mp1=lag.max+1
  
ACF1=stats::acf(y, lag.max, plot=FALSE)$acf[2:mp1]
PACF= stats::pacf(y, lag.max, plot=FALSE)$acf
  
LAG=1:lag.max
minA=min(ACF1);  minP=min(PACF)
maxA=max(ACF1);  maxP=max(PACF)
U=2/sqrt(num)
L=-U
# minu=min(minA,minP,L)-0.01
# minu = -1
minu=max(-1, min(minA,minP,L)-.01) ; maxu= min(1, max(maxA,maxP)+0.05)

op<-par(fig= c(0,1,0.6,1), mar=c(1,5,3,2))
plot(y, xlab="", cex=.8, cex.lab=0.9, ylab=names(y))
 par(new=TRUE, fig= c(0,1,0.32,0.55), mar=c(0,5,0,2))
  plot(LAG, ACF1, type="h", ylim=c(minu,maxu), xlab="", xaxt="n",
  cex=.8, cex.lab=0.9, las=1, ylab="ACF")
  abline(h=0)
  abline(h=L, lty="dashed", col="blue")
  abline(h=U, lty="dashed", col="blue")
  legend("topright", legend= paste("Time Series:",deparse(substitute(y))),
  bty="n", cex=.9)

  par(new=TRUE,fig= c(0,1,0,0.32), mar=c(3,5,0,2))
  plot(LAG, PACF, type="h" ,ylim=c(minu,maxu), xlab="", cex=.8, las=1,
  cex.lab=0.9, ylab="PACF")
  abline(h=0)
  abline(h=L, lty="dashed", col="blue")
  abline(h=U, lty="dashed", col="blue")
  if(numer)   return(cbind(LAG, ACF1, PACF))
par(op)  
}