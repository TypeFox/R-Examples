diffograph <-function (x, main=NULL,color1="red",color2="blue",color3="black",
cex.axis=0.8,las=1,pch=20,bty="l",cex=0.8,lwd=1,xlab="",ylab="",...)
{
  Ax<-x$means[,1]
  Ay<-x$comparison[,4:5]
  alpha<-x$parameters$alpha
  test<-x$parameters$test
  name.t<-x$parameters$name.t
  nombre<-rownames(x$means)
  n<-length(Ax)
  m<-nrow(Ay)
  k<-0
  p<-matrix(0,ncol=2,nrow=m)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      k<-k+1
      p[k,1]<-Ax[i]
      p[k,2]<-Ax[j]
    }
  }
  rango<-range(p)
  delta <- 0.3*abs(rango[2]-rango[1])
  if(rango[1]<0) rango[1]<-1.5*rango[1]
  if(rango[1]>0) rango[1]<-0.5*rango[1]
  if(rango[2]<0) rango[2]<-0.5*rango[2]
  if(rango[2]>0) rango[2]<-1.5*rango[2]
  if(rango[1]==0) rango[1]<- -delta
  if(rango[2]==0) rango[2]<-  delta
  d<-abs(Ay[,2]-Ay[,1])
  s<-Ay[,1]*Ay[,2]
  m<-length(p[,1])
  color<-rep(color1,m)
  for(k in 1:m){
    if(s[k]<=0) color[k]<-color2 
  }
  x1<-p[,1]-d/4
  y1<-p[,2]+d/4
  x2<-p[,1]+d/4
  y2<-p[,2]-d/4
  rxy<-range(c(x1,x2,y1,y2))
  par(mar=c(4,4,2,2))
  if(is.null(main))main=paste(names(x$means)[1],"Comparisons for",name.t)
  Az<-data.frame(name=nombre,value=Ax)
  Az<-Az[order(Az[,2]),]
  rownames(Az)<-Az[,2]
  impar<-seq(1,n,2)
  par<-seq(2,n,2)
  d<-matrix(c(1,2),ncol=1,nrow=2,byrow=T)
  nf<-graphics::layout(d, heights=c(6,1.5),respect=F)
  plot(p[,1],p[,2],xlim=rxy,ylim=rxy,cex.axis=cex.axis,las=las,cex=cex,
       xlab=xlab,ylab=ylab,bty=bty,main=main,cex.main=cex,...)
  segments(x1,y1,x2,y2,col=color,lwd=lwd)
  points(p[,1],p[,2],pch=pch,cex=cex,col=color)
  abline(a=0,b=1,col=color3,lty=2)
  text(Az[impar,2],rxy[1],Az[impar,1],cex=cex)
  text(Az[par,2],1.2*rxy[1],Az[par,1],cex=cex)
  text(rxy[2],Az[impar,2],Az[impar,1],cex=cex)
  text(0.96*rxy[2],Az[par,2],Az[par,1],cex=cex)
  abline(v=Az[,2],h=Az[,2],lty=4)
  par(mar=c(0,2,0,2))
  plot(1,1,axes=FALSE,cex=0)
  method<-paste('Differences for alpha =',alpha,"(",test,")")
  text(1,1,method,cex=cex)
    legend("bottom",c("Significant","Not significant"),lty=c(1,1),bty="n",
    lwd=2, cex=cex,  col=c(color1,color2),box.col ="white",horiz=TRUE)

}
