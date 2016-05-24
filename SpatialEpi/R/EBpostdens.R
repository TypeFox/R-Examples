EBpostdens <-
function(Y, E, alpha, beta, Xrow=NULL, lower=NULL, upper=NULL, main=""){
  if (is.null(Xrow)) Xrow <- matrix(c(1),nrow=1,ncol=1)
  xvals <- seq(lower,upper,.01)
  mu <- as.numeric(exp(Xrow %*% beta))
  yvals <- dgamma(xvals,alpha+Y,(alpha+E*mu)/mu)
  plot(xvals,yvals,type="n",
       xlab=expression(theta),ylab="EB density",cex.lab=1.2)
  title(paste(main))
  lines(xvals,yvals)
  lines(c(Y/E,Y/E),c(0,max(yvals)),lty=2)
}
