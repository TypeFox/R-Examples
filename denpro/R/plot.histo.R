plot.histo<-function(pcf,col=NULL,cex.axis=1,cex.lab=1,ylab="",xlab="",
xaxt="s",yaxt="s")
{
if (is.null(col)){
   f0<-sqrt(pcf$value)  #f0<-pcf$value
   colo<-1-(f0-min(f0)+0.5)/(max(f0)-min(f0)+0.5)
   #colo<-1-(f0-min(f0)+0.02)/(max(f0)+0.05-min(f0)+0.02)
   #colo<-1-(f0-min(f0))/(max(f0)-min(f0))
   col<-gray(colo)
}

d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]=(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i];

xmin<-pcf$support[1]
xmax<-pcf$support[2]
ymin<-pcf$support[3]
ymax<-pcf$support[4]

plot(xmin,ymin,type="n",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
pch=20,cex.axis=cex.axis,cex.lab=cex.lab,ylab=ylab,xlab=xlab,xaxt=xaxt,yaxt=yaxt)

lenni<-length(pcf$value)
for (i in 1:lenni){
     x1<-pcf$support[1]+step[1]*pcf$down[i,1]
     x2<-pcf$support[1]+step[1]*pcf$high[i,1] 
     y1<-pcf$support[3]+step[2]*pcf$down[i,2]
     y2<-pcf$support[3]+step[2]*pcf$high[i,2] 

     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=col[i],lty="blank")
}

}



