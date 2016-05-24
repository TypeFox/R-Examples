color.legend<-function (xl,yb,xr,yt,legend,rect.col,cex=1,align="lt",
 gradient="x",...) {

 oldcex<-par("cex")
 par(xpd=TRUE,cex=cex)
 gradient.rect(xl,yb,xr,yt,col=rect.col,nslices=length(rect.col),
  gradient=gradient)
 if(gradient == "x") {
  xsqueeze<-(xr-xl)/(2*length(rect.col))
  textx<-seq(xl+xsqueeze,xr-xsqueeze,length.out=length(legend))
  if(match(align,"rb",0)) {
   texty<-yb-0.2*strheight("O")
   textadj<-c(0.5,1)
  }
  else {
   # assume it's the default
   texty<-yt+0.2*strheight("O")
   textadj<-c(0.5,0)
  }
 }
 else {
  ysqueeze<-(yt-yb)/(2*length(rect.col))
  texty<-seq(yb+ysqueeze,yt-ysqueeze,length.out=length(legend))
  if(match(align,"rb",0)) {
   textx<-xr+0.2*strwidth("O")
   textadj<-c(0,0.5)
  }
  else {
   # assume it's the default
   textx<-xl-0.2*strwidth("O")
   textadj<-c(1,0.5)
  }
 }
 text(textx,texty,labels=legend,adj=textadj,...)
 par(xpd=FALSE,cex=oldcex)
}
