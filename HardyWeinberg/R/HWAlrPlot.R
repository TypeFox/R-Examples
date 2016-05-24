HWAlrPlot <- function(X,zeroadj=0.5) {
  # plots markers in ALR coordinates, with a line representing HWE.
  nullvalue <- (2/3)*log(2)
  Xalr <- HWAlr(X,zeroadj)
  opar <- par(mar=c(5,5,4,2)+0.1)
  plot(Xalr[,1],Xalr[,2],pch=19,xlab=expression(ln(f[AA]/f[AB])),ylab=expression(ln(f[BB]/f[AB])),asp=1,
       main="Alr-transformation",cex.main=2,cex.lab=2)
  vx <- c(par("usr")[1],par("usr")[2])
  vy <- c(nullvalue,nullvalue)
  p <- af(X)
  x <- log(p/(1-p))-log(2)
  y <- -log(p/(1-p))-log(2)
  lines(x,y)
  par(opar)
  return(NULL)
}
