tab.title<-function(label,text.col=par("fg"),tab.col=par("bg"),
 border=par("fg"),lwd=par("lwd"),cex=1.5,pad.mult=1.6,radius=0) {

 plotin<-par("pin")
 marin<-par("mai")
 plim<-par("usr")
 xmult<-(plim[2]-plim[1])/(plim[4]-plim[3])*plotin[2]/plotin[1]
 plim[3]<-plim[4]
 plim[4]<-plim[3]+((plim[2]-plim[1])/xmult)*marin[3]/plotin[2]
 oldcex<-par("cex")
 par(cex=cex)
 tabtop<-plim[3]+strheight(label)*pad.mult
 oldxlog<-par("xlog")
 oldylog<-par("ylog")
 par(xlog=FALSE,ylog=FALSE)
 if(radius == 0) {
  tabx<-c(plim[1],plim[1],plim[2],plim[2])
  taby<-c(plim[3],tabtop,tabtop,plim[3])
 }
 else {
  xradius<-radius*(tabtop-plim[3])*xmult
  xcurve1<-xradius*cos(seq(0,pi/2,length.out=20))
  xcurve2<-xradius*cos(seq(pi/2,pi,length.out=20))
  tabx<-c(plim[2],
   (plim[2]-xradius)+xcurve1,
   (plim[1]+xradius)+xcurve2,
   plim[1])
  yradius<-radius*(tabtop-plim[3])
  ycurve1<-yradius*sin(seq(0,pi/2,length.out=20))
  ycurve2<-yradius*sin(seq(pi/2,pi,length.out=20))
  taby<-c(plim[3],
   tabtop-yradius+ycurve1,
   tabtop-yradius+ycurve2,
   plim[3])
 }
 par(xpd=TRUE)
 polygon(tabx,taby,border=border,col=tab.col,lwd=lwd)
 text((plim[1]+plim[2])/2,(plim[3]+tabtop)/2,label,
  col=text.col)
 par(xpd=FALSE,xlog=oldxlog,ylog=oldylog,cex=oldcex)
}
