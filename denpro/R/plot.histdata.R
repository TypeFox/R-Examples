plot.histdata<-function(dendat,col,pcf,i1=1,i2=2,i3=3,
simple=FALSE,cut=dim(dendat)[1],
xlab="",xlim=NULL,ylim=NULL,cex.axis=1)
{
d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]=(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]

ord<-order(dendat[,i3])#,decreasing=TRUE)
ord<-ord[1:cut]
deudat<-dendat[ord,c(i1,i2)]

if (is.null(xlim)){
   xmin<-min(deudat[,1])
   xmax<-max(deudat[,1])
   xlim<-c(xmin,xmax)
}
if (is.null(ylim)){
   ymin<-min(deudat[,2])
   ymax<-max(deudat[,2])
   ylim=c(ymin,ymax)
}

if (simple) 
plot(deudat,pch=19,col=col[ord],xlab=xlab,ylab="",cex.axis=cex.axis,
xlim=xlim,ylim=ylim)

if (!simple){
pointx<-(xlim[1]+xlim[2])/2
pointy<-(ylim[1]+ylim[2])/2
plot(pointx,pointy,type="n",ylab="",xlim=xlim,ylim=ylim,,cex.axis=cex.axis,
pch=20,xlab=xlab)
for (i in 1:cut){
     mu1<-deudat[i,1]
     mu2<-deudat[i,2]
     x1<-mu1-step[i1]/2
     x2<-mu1+step[i1]/2
     y1<-mu2-step[i2]/2
     y2<-mu2+step[i2]/2
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=col[i])
}
}

}


