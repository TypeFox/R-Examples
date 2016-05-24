"edaplot" <-
function(data,scatter=TRUE,box=TRUE, P.plot=TRUE, D.plot=TRUE,
         P.main=paste("Histogram of",deparse(substitute(data))),
         P.sub=NULL, P.xlab=deparse(substitute(data)), P.ylab=default, P.ann=par("ann"),
         P.axes=TRUE, P.frame.plot=P.axes, P.log=FALSE, P.logfine=c(2,5,10), P.xlim=NULL,
         P.cex.lab=1.4,B.range=1.5, B.notch=FALSE, B.outline=TRUE,
         B.border=par("fg"), B.col=NULL, B.pch=par("pch"), B.cex=1, B.bg=NA, H.breaks="Sturges",
         H.freq=TRUE, H.include.lowest=TRUE, H.right=TRUE, H.density=NULL, H.angle=45,
         H.col=NULL, H.border=NULL, H.labels=FALSE, S.pch=".", S.col=par("col"), S.bg=NA, S.cex=1,
         D.lwd=1,D.lty=1)

{
# Combination of Histogram, density trace, 1-dim scattergram and boxplot

# P.log ... If TRUE: x-axis is in log-scale (Peter, 17.2.2006)
# P.logfine ... 10 for powers of 10, c(5,10) for powers of (5,10), c(2,5,10) for powers of (2,5,10)

 xhi<-hist(data,plot=FALSE, breaks=H.breaks, include.lowest=H.include.lowest, right=H.right)
 bxp<-boxplot(data,plot=FALSE, range=B.range, notch=B.notch, outline=B.outline)
 dens <- density(data)  # for density trace

 if ( P.plot )
 {
  x<-xhi$breaks
  if ( H.freq )
     {
      y<-xhi$counts
      default="Frequency"
      D.plot <- FALSE
     }
  else
     {
      y<-xhi$density
      default="Relative frequency"
     }
  h <- -max(y,dens$y)/8

  if ( scatter && box )
     {a<-2}
  else
  if ( (scatter & !box) | (!scatter && box))
     {a<-1}
  else
     {a<-0}

    if (D.plot) {   # plot density trace
      plot(x,c(y,a*h),type="n", main=P.main, sub=P.sub, xlab=P.xlab, ylab=P.ylab, ann=P.ann,
        frame.plot=P.frame.plot,ylim=c(min(c(y,a*h)),max(y,dens$y)),xaxt="n",yaxt="n",xlim=P.xlim,
        cex.lab=P.cex.lab)
      if (P.axes) {
        ay<-axTicks(2)
        axis(2,at=ay[ay>=0],labels=ay[ay>=0])
        if (P.log) axis(1,at=log10(alog<-sort(c((10^(-50:50))%*%t(P.logfine)))),labels=alog)
        else axis(1)
      } 
      lines(dens,lwd=D.lwd,lty=D.lty)
    }
    else {
      plot(x,c(y,a*h),type="n", main=P.main, sub=P.sub, xlab=P.xlab, ylab=P.ylab, ann=P.ann,
        axes=P.axes, frame.plot=P.frame.plot,xlim=P.xlim, cex.lab=P.cex.lab)
    }
  hist(data, add=TRUE, axes=FALSE, breaks=H.breaks, freq=H.freq, include.lowest=H.include.lowest,
      right=H.right, density=H.density, angle=H.angle, col=H.col, border=H.border, labels=H.labels)
    
  if ( scatter && box )
     {
      rect(xhi$breaks[1],h,xhi$breaks[length(xhi$breaks)],0)
      points(data,runif(length(data),0.1,0.9)*h,pch=S.pch, col=S.col, bg=S.bg, cex=S.cex)
      rect(xhi$breaks[1],2*h,xhi$breaks[length(xhi$breaks)],0)
      boxplot(data,add=TRUE,horizontal=TRUE,boxwex=-h*1.2,at=h*1.5, axes=FALSE, range=B.range, notch=B.notch,
              outline=B.outline, border=B.border, col=B.col, pch=B.pch, bg=B.bg, cex=B.cex,xlim=P.xlim)
     }

  else
  if ( scatter && !box )
     {
      rect(xhi$breaks[1],h,xhi$breaks[length(xhi$breaks)],0)
      points(data,runif(length(data),0.1,0.9)*h,pch=S.pch, col=S.col, bg=S.bg, cex=S.cex)
     }
  else
  if ( !scatter && box )
     {
      rect(xhi$breaks[1],h,xhi$breaks[length(xhi$breaks)],0)
      boxplot(data,add=TRUE,horizontal=TRUE,boxwex=-h*1.2,at=h*.5, axes=FALSE, range=B.range, notch=B.notch,
              outline=B.outline, border=B.border, col=B.col, pch=B.pch, bg=B.bg, cex=B.cex,xlim=P.xlim)
     }
 }

 l<-list(H=xhi, B=bxp)
 return(invisible(l))

}
