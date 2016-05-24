# To draw the sample paths
visualize.sp <- function(data,epsilon=0.05,nb.sp=10,plotfunc=plot,...) {
  
  limit <- max(abs(data[1:nb.sp,]))
  nmax <- ncol(data)
  plotfunc(data[1,],type="l",xlab="n",ylim=c(-limit,limit),ylab='',lwd=1,pch=46,...)
  mtext(bquote(bold("+" ~ epsilon)),side=2,line=2.5,xpd=TRUE,col="blue",at=epsilon,las=1)
  mtext(bquote(bold("-" ~ epsilon)),side=2,line=2.5,xpd=TRUE,col="blue",at=-epsilon,las=1)
  apply(data[1:nb.sp,],MARGIN=1,FUN=points,lwd=1,type="l",pch=46,...)
  abline(h=epsilon,col="blue")
  abline(h=-epsilon,col="blue")
  
  
}
