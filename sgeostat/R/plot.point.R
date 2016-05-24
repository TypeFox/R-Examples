# A plot function for the "point" class of objects
assign("plot.point",
function (x,v,legend.pos=0,axes=TRUE,xlab='',ylab='',add=FALSE,...) {
# Be careful to plot in a square region.  We can't distort the earth!
# I can't seem to force Splus to use the same scaling on both axes!
# But we must find a way!
  old.par <- par(pty='s')
  xdiff <- max(x$x) - min(x$x)
  ydiff <- max(x$y) - min(x$y)
  if (xdiff < ydiff) {
#   Set up our limits so that there are ydiff units on x and y...
    ylimits <- c(min(x$y),max(x$y))
    xlimits <- c((min(x$x) + xdiff/2) - ydiff/2,
                (min(x$x) + xdiff/2) + ydiff/2)
  }
  else {
    xlimits <- c(min(x$x),max(x$x))
    ylimits <- c((min(x$y) + ydiff/2) - xdiff/2,
                (min(x$y) + ydiff/2) + xdiff/2)
  }
  if (!missing(v)) {
    v.name<-v
    v <- x[[match(v,names(x))]]
#    colors <- cut(v,c(min(v)-1,quantile(v,c(.25,.5,.75,1))))
    v.q<-quantile(v,c(.25,.5,.75,1))	
    qcol<-c(3,7,2,6) # green, yellow, red, cyan
    colors <- cut(v,c(min(v)-1,v.q),labels=c(1:4))
    if(!add){
      plot(x$x,x$y, axes=axes,xlab=xlab,ylab=ylab,type='n',
           xlim=xlimits,ylim=ylimits)
    }
    for (i in as.numeric(unique(colors))){
      points(x$x[colors==i],x$y[colors==i],col=qcol[i],...)
}
    if (legend.pos!=0){
      l.x<-switch(legend.pos,
		  xlimits[1],xlimits[2],
		  xlimits[2],xlimits[1])
      l.xj<-switch(legend.pos,0,1,1,0)
      l.y<-switch(legend.pos,
		  ylimits[1],ylimits[1],
		  ylimits[2],ylimits[2])
      l.yj<-switch(legend.pos,0,0,1,1)
      legend(l.x,
	     l.y,
	     c(paste("[",min(v),",",v.q[1],"]"),
	       paste("(",v.q[1],",",v.q[2],"]"),
	       paste("(",v.q[2],",",v.q[3],"]"),
	       paste("(",v.q[3],",",max(v),"]")),
	     qcol,
	     xjust=l.xj,yjust=l.yj)
    }
    title(deparse(substitute(x)))
    
  }
  else {
    if(!add){
      plot(x$x,x$y, axes=axes,xlab=xlab,ylab=ylab,
           xlim=xlimits,ylim=ylimits)
      title(deparse(substitute(x)))
    } else {
      points(x$x,x$y)
    }
  }
  invisible(par(old.par))
})
