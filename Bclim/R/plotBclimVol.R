plotBclimVol <-
function(x,dim=1,title=NULL,presentleft=TRUE,denscol="red",denstransp=0.5,leg=TRUE,mean=TRUE,legloc="topleft",...) {
  
  if(class(x)!="Bclim") stop("Needs a Bclim output object")
  
  # Create HDRs for volatility
  errorbar <- matrix(NA,nrow=length(x$time.grid)-1,ncol=3)
  for(i in 1:(length(x$time.grid)-1)) errorbar[i,] <- quantile(x$v.interp[,i,dim],probs=c(0.025,0.5,0.975))
  for(i in 1:(length(x$time.grid)-1)) errorbar[i,2] <- mean(x$v.interp[,i,dim])
  
  # Sort out colours
  tmp <- col2rgb(denscol)
  mycol <- rgb(tmp[1,1]/255,tmp[2,1]/255,tmp[3,1]/255)
  mycol2 <- paste(mycol,as.character(as.hexmode(round(denstransp*255,0))),sep="")
  
  # Set up plot 
  par(mar=c(4,4,3,1))
  xrange <- range(c(0,x$time.grid))
  if(!presentleft) xrange <- rev(xrange)
  yrange <- range(c(0,as.vector(errorbar)))
  mytitle <- title
  if(is.null(title)) mytitle <- paste(x$core.name,": ",x$clim.dims[dim],sep="")
  
  if(dim==1) {
    plot(1,1,type="n",xlim=xrange,ylim=yrange,xaxt='n',xlab="Age (k cal years BP)",ylab=expression(paste("GDD5 volatility (",degree,"C days)",sep="")),las=1,bty="n",main=mytitle)
  }
  if(dim==2) {
    plot(1,1,type="n",xlim=xrange,ylim=yrange,xaxt='n',xlab="Age (k cal years BP)",ylab=expression(paste("MTCO volatility (",degree,"C)",sep="")),las=1,bty="n",main=mytitle)
  }
  if(dim==3) {
    plot(1,1,type="n",xlim=xrange,ylim=yrange,xaxt='n',xlab="Age (k cal years BP)",ylab=x$clim.dims[dim],las=1,bty="n",main=mytitle)
  }
  
  axis(side=1,at=pretty(x$time.grid,n=10))
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="lightgray",border="NA")
  grid(col="white")  
  
  # Draw lines    
  tgridshift = x$time.grid[-1]-c(diff(x$time.grid)/2)
  lines(c(min(x$time.grid),max(x$time.grid)),c(min(yrange),min(yrange)))
  for(i in 1:length(x$time.grid[-1])) {
    lines(c(tgridshift[i],tgridshift[i]),c(errorbar[i,1],errorbar[i,3]),col=denscol)
    if(mean) points(tgridshift[i],errorbar[i,2],col=denscol)
  }
  for(i in 1:length(x$time.grid)) lines(c(x$time.grid[i],x$time.grid[i]),c(min(yrange),0.05*max(yrange)))
  
  # Finally draw a legend
  if(leg==TRUE) {
    if(!mean) legend(legloc,legend=c("95% credibility interval",'Time Interval'),lty=c(1,1),col=c(mycol2,'black'),pch=c(-1,-1),bty="n")
    if(mean) legend(legloc,legend=c("95% credibility interval","Mean",'Time Interval'),lty=c(1,-1,1),col=c(mycol2,denscol,'black'),pch=c(-1,1,-1),bty="n")
    #fill=c(denscol,NULL)
  }
    
# End of function   
}
