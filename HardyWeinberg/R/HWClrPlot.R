HWClrPlot <- function(X,zeroadj=0.5) {
  # makes a plot of the rows of X in centred log-ratio coordinates, and adds a line representing HWE.
  nullvalue <- (2/3)*log(2)
  Xclr <- HWClr(X,zeroadj)
  opar <- par(mar=c(5,5,4,2)+0.1)
  plot(Xclr[,1],Xclr[,2],pch=19,xlab=expression(ln(f[AA]/g)),ylab=expression(ln(f[AB]/g)),asp=1,main="Clr-transformation",cex.main=2,cex.lab=2)
  vx <- c(par("usr")[1],par("usr")[2])
  vy <- c(nullvalue,nullvalue)
  lines(vx,vy)
  par(opar)
  return(NULL)
}
