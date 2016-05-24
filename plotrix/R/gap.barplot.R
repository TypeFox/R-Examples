gap.barplot<-function (y,gap,xaxlab,xtics,yaxlab,ytics,xlim=NA,ylim=NA,
 xlab=NULL,ylab=NULL,halfwidth=NA,horiz=FALSE,col=NULL,...) {
 if (missing(y)) stop("y values required")
 if(missing(xtics)) xtics <- 1:length(y)
 if (missing(gap)) stop("gap must be specified")
 if (is.null(ylab)) ylab <- deparse(substitute(y))
 if (is.null(col)) col <- color.gradient(c(0,1),c(0,1,0),c(1,0),length(y))
 else if(length(col) < length(y)) rep(col,length.out=length(y))
 littleones <- which(y <= gap[1])
 bigones <- which(y >= gap[2])
 valid.y<-y[!is.na(y)]
 if(any(valid.y > gap[1] & valid.y < gap[2]))
  warning("gap includes some values of y")
 gapsize <- gap[2] - gap[1]
 if(missing(xaxlab)) xaxlab <- as.character(xtics)
 if(is.na(xlim[1])) xlim <- range(xtics)
 if(is.na(ylim[1])) ylim <- c(min(valid.y),max(valid.y) - gapsize)
 if(missing(ytics)) ytics <- pretty(y)
 if(missing(yaxlab)) yaxlab <- ytics
 littletics <- which(ytics < gap[1])
 bigtics <- which(ytics >= gap[2])
 if(is.na(halfwidth)) halfwidth <- min(diff(xtics))/2
 if(horiz) {
  if(!is.null(xlab)) {
   tmplab<-xlab
   xlab<-ylab
   ylab<-tmplab
  }
  plot(0,xlim=ylim,ylim=xlim,xlab=xlab,ylab=ylab,axes=FALSE,type="n",...)
  plot.lim <- par("usr")
  botgap<-ifelse(gap[1]<0,gap[1],xlim[1])
  box()
  axis(2,at=xtics,labels=xaxlab,...)
  axis(1,at=c(ytics[littletics],ytics[bigtics]-gapsize),
   labels=c(yaxlab[littletics],yaxlab[bigtics]),...)
  rect(botgap,xtics[y<gap[1]] - halfwidth,y[y<gap[1]],
   xtics[y<gap[1]] + halfwidth,col=col[y<gap[1]])
  rect(botgap,xtics[bigones] - halfwidth,y[bigones]-gapsize,
   xtics[bigones] + halfwidth,col=col[bigones])
  axis.break(1,gap[1],style="gap")
 }
 else {
  plot(0,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,axes=FALSE,type="n",...)
  plot.lim <- par("usr")
  botgap<-ifelse(gap[1]<0,gap[1],ylim[1])
  box()
  axis(1,at=xtics,labels=xaxlab,...)
  axis(2,at=c(ytics[littletics],ytics[bigtics] - gapsize),
   labels=c(yaxlab[littletics],yaxlab[bigtics]),...)
  rect(xtics[y<gap[1]] - halfwidth,botgap,xtics[y<gap[1]] + halfwidth,
   y[y<gap[1]],col=col[y<gap[1]])
  rect(xtics[bigones] - halfwidth,botgap,xtics[bigones] + halfwidth,
   y[bigones]-gapsize,col=col[bigones])
  axis.break(2,gap[1],style="gap")
 }
 invisible(xtics)
}
