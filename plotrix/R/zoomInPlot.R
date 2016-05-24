zoomInPlot<-function(x,y=NULL,xlim=NULL,ylim=NULL,rxlim=xlim,rylim=ylim,
 xend=NA,zoomtitle=NULL,titlepos=NA,...) {
 par(mfrow=c(1,2))
 if(is.null(y)) {
  y<-x
  x<-1:length(x)
 }
 if(is.null(xlim[1])) xlim<-range(x)
 if(is.null(ylim[1])) ylim<-range(y)
 plot(x,y,xlim=xlim,ylim=ylim,...)
 xext<-yext<-rxext<-ryext<-0
 if(par("xaxs") == "r") {
  xext<-diff(xlim)*0.04
  rxext<-diff(rxlim)*0.04
 }
 if(par("yaxs") == "r") {
  yext<-diff(ylim)*0.04
  ryext<-diff(rylim)*0.04
 }
 if(is.na(rxlim[1])) {
  newbox<-locator(2)
  rxlim<-sort(newbox$x)
  rylim<-sort(newbox$y)
 }
 rect(rxlim[1]-rxext,rylim[1]-ryext,rxlim[2]+rxext,rylim[2]+ryext)
 xylim<-par("usr")
 xypin<-par("pin")
 rxi0<-xypin[1]*(xylim[2]-(rxlim[1]-rxext))/diff(xylim[1:2])
 rxi1<-xypin[1]*(xylim[2]-(rxlim[2]+rxext))/diff(xylim[1:2])
 y01i<-xypin[2]*(xylim[4]-(rylim[2]+ryext))/diff(xylim[3:4])
 y02i<-xypin[2]*((rylim[1]-ryext)-xylim[3])/diff(xylim[3:4])
 plot(x,y,xlim=rxlim,ylim=rylim,...)
 xypin<-par("pin")
 par(xpd=NA)
 xylim<-par("usr")
 xymai<-par("mai")
 x0<-xylim[1]-diff(xylim[1:2])*(xymai[2]+xymai[4]+rxi0)/xypin[1]
 x1<-xylim[1]-diff(xylim[1:2])*(xymai[2]+xymai[4]+rxi1)/xypin[1]
 y01<-xylim[4]-diff(xylim[3:4])*y01i/xypin[2]
 y02<-xylim[3]+diff(xylim[3:4])*y02i/xypin[2]
 par(xpd=TRUE)
 if(!is.null(zoomtitle)) {
  if(is.na(titlepos)) titlepos<-getFigCtr()[1]
  mtext(zoomtitle,at=titlepos,cex=1.5,line=1.5)
 }
 if(is.na(xend)) xend<-xylim[1]-diff(xylim[1:2])*xymai[2]/(2*xypin[1])
 xprop0<-(xylim[1]-xend)/(xylim[1]-x0)
 xprop1<-(xylim[2]-xend)/(xylim[2]-x1)
 par(xpd=NA)
 segments(c(x0,x0,x1,x1),c(y01,y02,y01,y02),
  c(xend,xend,xend,xend),
  c(xylim[4]-(xylim[4]-y01)*xprop0,xylim[3]+(y02-xylim[3])*xprop0,
  xylim[4]-(xylim[4]-y01)*xprop1,xylim[3]+(y02-xylim[3])*xprop1))
 par(xpd=FALSE,mfrow=c(1,1))
}
