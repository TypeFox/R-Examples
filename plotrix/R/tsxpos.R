tsxpos<-function(x,xlim,nint) {
 # make sure that there is a graphics device open
 if(dev.cur() == 1) stop("A graphics device must be open")
 if(missing(xlim)) {
  plotlim<-par("usr")
  divlim<-ifelse(par("xaxs") == "r",0.04,0)
  xrange<-(plotlim[2]-plotlim[1])/(1+2*divlim)
  xstart<-plotlim[1]+xrange*divlim
 }
 else {
  xstart<-xlim[1]
  xrange<-diff(xlim)
 }
 if(missing(nint)) nint<-length(x)-1
 return(cumsum(c(xstart,rep(xrange/nint,nint))))
}
