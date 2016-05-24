dhist <- function(x,fac,col,legend,pos.legend,title.legend=NULL,lab.legend=NULL,xlab,ylab=NULL,
  drawextaxes=TRUE,drawintaxes=FALSE,xlim=NULL,...) {
  ymax <- integer(nlevels(fac))
  for (i in 1:nlevels(fac)) {
    ymax[i] <- max(density(x[as.numeric(fac)==i])$y)
  }
  h <- suppressWarnings(hist(x,freq=FALSE,plot=FALSE))
  oldmar <- par()$mar
  if (is.null(ylab)) {ylab="Probability density"}
  if (!drawextaxes) {par(mar=c(3.1,2.1,2.1,0.1))}
  xlim <- if(!is.null(xlim)) {xlim} else {range(h$breaks)}
  plot(0,xlim=xlim,ylim=c(0,max(ymax)),xlab="",ylab="",cex=0,axes=FALSE,...)
  if(drawextaxes) {
    axis(1)
    axis(2)
  }
  if (drawintaxes) {abline(v=0,col="grey")}
  lab.line <- c(ifelse(drawextaxes,3,1.2),ifelse(drawextaxes,3,0.6))
  mtext(c(xlab,ylab),side=c(1,2),line=lab.line,at=c(mean(range(x)),mean(c(0,max(ymax)))))
  dens <- tapply(x,fac,function(x) density(x))
  if (!is.numeric(col)) {
    col3 <- col4 <- col
  } else {
    col2 <- col2rgb(palette()[col])
    col3 <- apply(col2,2,function(x) rgb(x[1],x[2],x[3],alpha=0.4*255,maxColorValue=255))
    col4 <- apply(col2,2,function(x) rgb(x[1],x[2],x[3],alpha=255,maxColorValue=255))  
  }
  for (i in 1:nlevels(fac)) {
    d <- dens[[i]]
    polygon(d$x,d$y,col=col3[i],border=NA)
    rug(x[as.numeric(fac)==i],col=col4[i])
  }
  box()
  if (legend) {
    if (is.null(lab.legend)) {lab.legend <- levels(fac)}
    if (!is.null(title.legend) && nchar(title.legend)>0) {
	legend(pos.legend,lab.legend,fill=col3,title=title.legend)
    } else {
	legend(pos.legend,lab.legend,fill=col3)
    }
  }
  par(mar=oldmar)
}
