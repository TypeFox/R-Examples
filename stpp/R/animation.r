animation <- function(xyt, s.region, t.region, runtime=1, incident="red", prevalent="pink3", pch=19, cex=0.5, plot.s.region=TRUE, scales=TRUE, border.frac=0.05, add=FALSE)
{ 
  if (missing(s.region)) s.region <- sbox(xyt[,1:2],xfrac=0.01,yfrac=0.01)
  if (missing(t.region)) t.region <- range(xyt[,3],na.rm=T)
  
  ott<-order(xyt[,3])
  sxyt<-xyt[ott,]
  rangex<-range(s.region[,1])
  rangey<-range(s.region[,2])
  xlim<-c(rangex[1]-border.frac*(rangex[2]-rangex[1]),rangex[2]+border.frac*(rangex[2]-rangex[1]))
  ylim<-c(rangey[1]-border.frac*(rangey[2]-rangey[1]),rangey[2]+border.frac*(rangey[2]-rangey[1]))
  xy<-as.matrix(sxyt[,1:2])
  tt<-sxyt[,3]
  npts<-length(tt)
  T0 <- max(t.region)

  if (add==FALSE)
    {
	par(pty="s",mfrow=c(1,1))
      if (scales==FALSE)
        plot(xy[,1],xy[,2],type="n",xlim=xlim,ylim=ylim,xaxt="n",yaxt="n",bty="n",xlab=" ",ylab=" ")
      if (scales==TRUE)
        plot(sxyt[,1],sxyt[,2],type="n",xlim=xlim,ylim=ylim,bty="n",xlab="X",ylab="Y")
      if (plot.s.region==TRUE)
        polymap(as.points(s.region),add=TRUE,lwd=2)
    }
  nplotted<-0
  tt.now<-0
  tt.fade<-0	
  while (nplotted<npts)
    {
      i<-nplotted+1
      tt.gap<-tt[i]-tt.now
	if (!(is.null(runtime)))
      	junk<-Sys.sleep((tt.gap/T0)*runtime)
      n.fade<-sum(tt[1:i]<=tt.fade)
      if (sum(n.fade)>0)
        points(xy[1:n.fade,1],xy[1:n.fade,2],col=prevalent,pch=19,cex=cex)

	M=which(tt==tt[i])
	points(xy[M,1],xy[M,2],col=incident,pch=19,cex=cex)
      nplotted<-nplotted+length(M)
      tt.now<-tt[i]
      tt.fade<-tt.now
    }
}


