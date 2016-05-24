plot.cv.sparsenet=function(x,...){
  oldpar=par(mar=c(4,4,3,5))
  on.exit(par(oldpar))
  parms=x$sparsenet.fit$parms
  cvm=x$cvm
  ngamma=ncol(cvm)
  cvcolors=rainbow(ngamma,start=.1,end=1)
  plot(log(t(parms[2,,])),x$cvm,type="n",xlab=expression(log(lambda)),ylab="CV Error")
  for( i in seq(ngamma)) lines(log(parms[2,i,]),x$cvm[,i],lwd=2,col=cvcolors[i])
  mins=c(log(x$parms.min[2]),log(x$parms.1se[2]))
  abline(v=mins,lty=2)
  nzero=x$nzero
  ccds=rbind(x$which.min,x$which.1se)
  dof=nzero[ccds]
  axis(3,at=mins,labels=format(dof),line=0)
  
  colorlegend(posy=c(0.2,.8),posx=c(.93,.945),col=rev(rainbow(ngamma,start=.1,end=1)),zlim=c(0,ngamma),zval=c(0,1,3,5,ngamma),main=expression(log ( gamma )))
  invisible()
}
