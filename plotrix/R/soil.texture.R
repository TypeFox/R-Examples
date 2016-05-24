get.soil.texture<-function(use.percentages=FALSE,cnames=c("sand","silt","clay")) {
 return(get.triprop(use.percentages=use.percentages, cnames=cnames))
}

soil.texture<-function(soiltexture=NULL,main="",at=seq(0.1,0.9,by=0.1),
 axis.labels=c("percent sand","percent silt","percent clay"),
 tick.labels=list(l=seq(10,90,by=10),r=seq(10,90,by=10),b=seq(10,90,by=10)),
 show.names=TRUE,show.lines=TRUE,col.names="gray",bg.names=par("bg"),
 show.grid=FALSE,col.axis="black",col.lines="gray",col.grid="gray",
 lty.grid=3,show.legend=FALSE,label.points=FALSE,point.labels=NULL,
 col.symbols="black",pch=par("pch"), ...) {

 par(xpd=TRUE)
 plot(0.5,type="n",axes=FALSE,xlim=c(0,1),ylim=c(0,1),main=NA,xlab=NA, ylab=NA)
 triax.plot(x=NULL,main=main,at=at,axis.labels=axis.labels,tick.labels=tick.labels,
  col.axis=col.axis,show.grid=show.grid,col.grid=col.grid,lty.grid=lty.grid)
 arrows(0.12,0.41,0.22,0.57,length=0.15)
 arrows(0.78,0.57,0.88,0.41,length=0.15)
 arrows(0.6,-0.1,0.38,-0.1,length=0.15)
 if(show.lines) {
  triax.segments<-function(h1,h3,t1,t3,col) {
   segments(1-h1-h3/2,h3*sin(pi/3),1-t1-t3/2,t3*sin(pi/3),col=col)
  }
  # from bottom-left to up
  h1 <- c(85, 70, 80, 52, 52, 50, 20,  8, 52, 45, 45, 65, 45, 20, 20)/100
  h3 <- c( 0,  0, 20, 20,  7,  0,  0, 12, 20, 27, 27, 35, 40, 27, 40)/100
  t1 <- c(90, 85, 52, 52, 43, 23,  8,  0, 45,  0, 45, 45,  0, 20,  0)/100
  t3 <- c(10, 15, 20,  7,  7, 27, 12, 12, 27, 27, 55, 35, 40, 40, 60)/100
  triax.segments(h1, h3, t1, t3, col.lines)
 }
 if(show.names) {
  xpos <- c(0.5, 0.7, 0.7, 0.73, 0.73, 0.5, 0.275, 0.275, 0.27,
   0.27, 0.25, 0.135, 0.18, 0.055, 0.49, 0.72, 0.9)
  ypos <- c(0.66, 0.49, 0.44, 0.36, 0.32, 0.35, 0.43, 0.39, 0.3,
   0.26, 0.13, 0.072, 0.032, 0.024, 0.18, 0.15, 0.06)*sin(pi/3)
  snames <- c("clay", "silty", "clay", "silty clay", "loam",
   "clay loam", "sandy", "clay", "sandy clay", "loam", "sandy loam",
   "loamy", "sand", "sand", "loam", "silt loam", "silt")
  boxed.labels(xpos, ypos, snames, border=FALSE, xpad=0.5)
 }
 par(xpd=FALSE)
 # now call triax.points
 if(is.null(soiltexture)) return(NULL)
 soilpoints <- triax.points(soiltexture, show.legend=show.legend,
  label.points=label.points, point.labels=point.labels,
  col.symbols=col.symbols, pch=pch, ...)
 invisible(soilpoints)
}
