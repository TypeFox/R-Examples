geometric_distance <- function(x,y,cx,cy,semi.major,semi.minor,rote.rad,ti,pred.method,ampx,ampy,lag,period) {
#if (pred.method=="find.times") {
  ti<-numeric(length(x))
for (i in 1:length(x)) {
  x0<-x[i]
  y0<-y[i]
  zmin1<-optimize(ellipsespot,c(0,pi),"x0"=x0,"y0"=y0,"cx"=cx,"cy"=cy,"semi.major"=semi.major,"semi.minor"=semi.minor,"rote.rad"=rote.rad)
  zmin2<-optimize(ellipsespot,c(pi,2*pi),"x0"=x0,"y0"=y0,"cx"=cx,"cy"=cy,"semi.major"=semi.major,"semi.minor"=semi.minor,"rote.rad"=rote.rad)
  ti[i]<-ifelse(zmin1$objective < zmin2$objective, zmin1, zmin2)[[1]]
  }
  pred.x<-cx +semi.major*cos(rote.rad)*cos(ti)-semi.minor*sin(rote.rad)*sin(ti)
  pred.y<-cy +semi.major*sin(rote.rad)*cos(ti)+semi.minor*cos(rote.rad)*sin(ti)
#}
#else {
#  srad<-sin(ti)
#  crad <- cos(ti)
#  ZZ <- cbind(srad,crad,rep(1,length(srad)))
#  phase.Ta.nls.polar <- atan2((- c(1,0,0)%*% solve(t(ZZ)%*%ZZ)%*%t(ZZ)%*%x),(c(0,1,0)%*% solve(t(ZZ)%*%ZZ)%*%t(ZZ)%*%x))
#  pred.x <- ampx*cos(ti+phase.Ta.nls.polar)+cx
#  pred.y <- ampy*cos(ti+phase.Ta.nls.polar-lag*2*pi/period)+cy
#}

return(list("pred.x"=pred.x,"pred.y"=pred.y,"period.time"=ti))
}
