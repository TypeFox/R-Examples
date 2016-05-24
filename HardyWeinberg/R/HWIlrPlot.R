HWIlrPlot <- function(X,zeroadj=0.5,...) {
  # Plot the rows of X in isometric logratio coordinates.
  nullvalue <- -sqrt(2/3)*log(2)
  Xilr <- HWIlr(X,zeroadj)
  opar <- par(mar=c(5,5,2,1))
  plot(Xilr[,1],Xilr[,2],pch=19,xlab=expression((1/sqrt(2))*ln(f[AA]/f[BB])),
       ylab=expression((1/sqrt(6))*ln(f[AA]*f[BB]/f[AB]^2)),...)
  vx <- c(par("usr")[1],par("usr")[2])
  vy <- c(nullvalue,nullvalue)
  lines(vx,vy)
  par(opar)
}
