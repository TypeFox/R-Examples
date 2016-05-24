plot.fitmarkov <-
function(x,...) {
     zz<- data.frame(x$fitted.data)
     yy<- data.frame(x$raw.data)
     p <- length(yy[1,]) 
#     check.ltypes<- is.null(l.types)
#     if(check.ltypes == TRUE) {
     l.types<- rep(c(1,1,1,1,1,1,1,3,3,3,3,3,3,3),ceiling(p/14))
#     }
     t<- x$t.measured
     s.types<- seq(1,p,1)
     tmod<- x$t.modeled
#
# line widths, types and colors
#
# The next 2 lines set defaults
   l.widths<- rep(c(3.6,2.8,2.3,1.9,1.5,1.1,0.6),ceiling(p/7))
   colors<- rgb(0,0,0,rep(c(35,70,100,135,175,200,240),ceiling(p/7)),maxColorValue=255)
# symbol types s.types
   s.types<- seq(1,p,1)
#
# raw plot
# Upper graph
  plot(c(min(t),max(t)),c(0,1.0),type="n",xlab="Year",ylab="Relative cover",cex.lab=0.8,main="",cex.axis=0.8)
  for (i in 1:p) lines(t,yy[,i],col=colors[i],lty=l.types[i],lwd=l.widths[i])                 # plot raw data
  for (i in 1:p) points(t,yy[,i],col=colors[i],pch=s.types[i],cex=0.6)                        # plot raw data
  legend("top",names(yy[1:p]),lty=l.types,pch=s.types,ncol=3,bty="n",pt.lwd=1,lwd=l.widths,cex=0.7,col=colors)
# Lower graph
  plot(c(min(t),max(t)),c(0,1.0),type="n",xlab="Year",ylab="Relative cover",cex.lab=0.8,main="",cex.axis=0.8)
  for (i in 1:p) lines(tmod,zz[,i],col=colors[i],lty=l.types[i],lwd=l.widths[i])              # plot simulated data
  for (i in 1:p) points(tmod,zz[,i],col=colors[i],pch=s.types[i],cex=0.6)              # plot simulated data
  legend("top",names(zz[1:p]),lty=l.types,pch=s.types,ncol=3,bty="n",pt.lwd=1,lwd=l.widths,cex=0.7,col=colors)
}
