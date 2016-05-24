checking.plots <-
function(model,n.id=3,COL=c("#0080FF","#A9E2FF"))
{
  Identify <- function(x,y,n.id=3,AD=0.2)
    {
    y[is.na(y)]<-0
    n <- length(y)
    y[is.infinite(y)]<-0
    pos <-as.numeric(names(y))
    oy <- order(abs(y))
    b <- y*AD
    W <- oy[(n - n.id + 1):n]
    text(x[W], y[W]+b[W], as.character(pos[W]))
    }
  par(mfrow=c(1,3),pty="s")
  require(MASS)
  varmodel <- deparse(substitute(model))
  y <- stdres(model)
  mAx <- max(abs(y[is.finite(y)])) + .5
  n <- length(y)
  plot( (1:n), y, ylab="standardized residuals",col=COL[1],
  xlab="ordered values",ylim=c(-mAx,mAx),
  main=paste("Standardized residuals versus \n ordered values for",varmodel))
  Identify((1:n),y,n.id)
  abline(h=c(2,3),col=COL,lty=2)
  abline(h=c(-2,-3),col=COL,lty=2)
  qqnorm(y,col=COL[1],ylab="standardized residuals",xlim=c(-mAx,mAx),
  ylim=c(-mAx,mAx),main=paste("Normal Q-Q plot of standardized \n residuals from",varmodel))
  abline(a=0,b=1,col=COL[2])
  junk <- qqnorm(y,plot.it=FALSE)
  Identify(junk$x,junk$y,n.id)
  plot(fitted(model),y,col=COL[1],xlab="fitted values",ylab="standardized residuals",
  ylim=c(-mAx,mAx), main=paste("Standardized residuals versus \n fitted values for",varmodel))
  abline(h=0,lty=2,col=COL[2])
  abline(h=c(2,3),col=COL,lty=2)
  abline(h=c(-2,-3),col=COL,lty=2)
  Identify(fitted(model),y,n.id)
  par(mfrow=c(1,1),pty="m")
}

