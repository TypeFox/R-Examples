magcon <-
function(x,y,h,doim=TRUE,docon=TRUE,dobar=TRUE,n=100,add=FALSE,xlab='',ylab='',imcol=rev(rainbow(1000,start=0,end=2/3)),conlevels=c(0.5,pnorm(1)-pnorm(-1),0.95), barposition='topright', barorient='v',bartitle='Contained %',bartitleshift=0,xlim=NULL,ylim=NULL,weights=NA,...){
if(is.null(xlim)){xlim=range(x,na.rm=TRUE)}
if(is.null(ylim)){ylim=range(y,na.rm=TRUE)}
use=x>=min(xlim) & x<=max(xlim) & y>=min(ylim) & y<=max(ylim)
if(is.na(weights[1])==FALSE & length(weights)==length(x)){weights=weights[use]}
x=x[use];y=y[use]
conlevels=1-conlevels
tempcon=sm.density(cbind(x,y),h=h,weights=weights,display='none',ngrid=n)
tempcon$x=tempcon$eval.points[,1]
tempcon$y=tempcon$eval.points[,2]
tempcon$z=tempcon$estimate
temp=sort(tempcon$z)
tempsum=cumsum(temp)
convfunc=approxfun(tempsum,temp)
levelmap=approxfun(convfunc(seq(0,1,len=1000)*max(tempsum)),seq(0,1,len=1000))
tempcon$z=matrix(levelmap(tempcon$z),nrow=n)
tempcon$z[is.na(tempcon$z)]=min(tempcon$z,na.rm=TRUE)
if(doim){
  if(add==FALSE){
    plot.new()
    plot.window(xlim=xlim,ylim=ylim)
    usrlims=par()$usr
    rect(usrlims[1],usrlims[3],usrlims[2],usrlims[4],col=imcol[1])
  }
  image(tempcon,col=imcol,axes=FALSE,add=TRUE,xlim=xlim,ylim=ylim)
}
if(doim & docon){contour(tempcon,levels=conlevels,add=TRUE,drawlabels=F,axes=FALSE,...)}
if(doim==FALSE & docon){contour(tempcon,levels=conlevels,add=add,drawlabels=F,axes=FALSE,xlim=xlim,ylim=ylim,...)}
if(add==FALSE){magaxis(xlab=xlab,ylab=ylab)}
if(doim & dobar){magbar(position=barposition,range=c(0,100),orient=barorient,col=rev(imcol),title=bartitle,titleshift=bartitleshift)}
return=tempcon
}
