draw.levset<-function(pcf,lev=NULL,bary=NULL,propor=0.1,col=NULL,
bound=NULL,dendat=NULL,xaxt="s",yaxt="s",cex.axis=1)
{

if (is.null(lev)) lev<-propor*max(pcf$value)

d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]=(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i];

if (is.null(bound)){
  xmin<-pcf$support[1]
  xmax<-pcf$support[2]
  ymin<-pcf$support[3]
  ymax<-pcf$support[4]
}
else{
  xmin<-bound[1]
  xmax<-bound[2]
  ymin<-bound[3]
  ymax<-bound[4]
}

if (is.null(bary))
   plot(xmin,ymin,type="n",xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
   pch=20,col="red",xaxt=xaxt,yaxt=yaxt,cex.axis=cex.axis)
else
  plot(x=bary[1],y=bary[2],
  xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
  pch=20,col="red",xaxt=xaxt,yaxt=yaxt,cex.axis=cex.axis)

lenni<-length(pcf$value)
for (i in 1:lenni){
  if (pcf$value[i]>=lev){

     x1<-pcf$support[1]+step[1]*pcf$down[i,1]
     x2<-pcf$support[1]+step[1]*pcf$high[i,1] 
     y1<-pcf$support[3]+step[2]*pcf$down[i,2]
     y2<-pcf$support[3]+step[2]*pcf$high[i,2] 

     if (is.null(col)) polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2))  
     else  polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=col[i],lty="blank")
  }
  i<-i+1
}

if (!is.null(dendat)) points(dendat)

#points(x=bary[1],y=bary[2],pch=20,col="red")

}



