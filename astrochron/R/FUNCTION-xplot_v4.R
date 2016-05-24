### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### xplot: create cross plot of two variable, with density estimates along
###        axes - (SRM: April 16, 2014; February 7, 2015; February 15, 2015;
###                     March 19, 2015)
###
###########################################################################


xplot <- function (x,y,xlab=NULL,ylab=NULL,main=NULL,fill=T)
{

  cat("\n----- GENERATING CROSS-PLOT WITH DENSITY ESTIMATES ON AXES -----\n")
  x<-data.frame(x)
  y<-data.frame(y)
  if( length(x[,1]) != length(y[,1]) ) stop("***** ERROR: number of data points in two series is not equivalent!")

  xmax=max(x[,1])
  xmin=min(x[,1])
  ymax=max(y[,1])
  ymin=min(y[,1])

  dev.new(title=paste("xplot results"),height=6,width=5.6)
  par(fig=c(0,0.8,0,0.8))
  plot(x[,1],y[,1],xlim=c(xmin,xmax),ylim=c(ymin,ymax),cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="",ylab="")
  if(is.null(xlab)) xlab=colnames(x)
  if(is.null(ylab)) ylab=colnames(y)
  if(is.null(main)) main=""
  mtext(ylab,side=2,line=2.5,cex=1.2)
  mtext(xlab,side=1,line=2.5,cex=1.2)
  par(fig=c(0,0.8,0.55,1), new=TRUE)
  z <- density(x[,1], bw="nrd0")
  plot(z$x,z$y, type='l',axes=FALSE,xlab="",ylab="",lwd=3,xlim=c(xmin,xmax))
  if(fill) 
   {
# find minimum (for plotting polygon)
    minY=min(z$y)
    polyY=append(minY,z$y)
    polyY=append(polyY,minY)
    minX=min(z$x)
    minX=minX-(minX*10^-10)
    maxX=max(z$x)
    maxX=maxX+(maxX*10^-10)
    polyX=append(minX,z$x)
    polyX=append(polyX,maxX)  
    polygon(polyX,polyY,col="#BEBEBE5A",border=NA)
   }
  par(fig=c(0.65,1,0,0.8),new=TRUE)
  z <- density(y[,1], bw="nrd0")
  plot(z$y, z$x, type='l',axes=FALSE,xlab="",ylab="",lwd=3,ylim=c(ymin,ymax))
  if(fill) 
   {
# find minimum (for plotting polygon)
    minY=min(z$y)
    polyY=append(minY,z$y)
    polyY=append(polyY,minY)
    minX=min(z$x)
    minX=minX-(minX*10^-10)
    maxX=max(z$x)
    maxX=maxX+(maxX*10^-10)
    polyX=append(minX,z$x)
    polyX=append(polyX,maxX)  
    polygon(polyY,polyX,col="#BEBEBE5A",border=NA)
   }
  mtext(main, side=3, cex=1.5, outer=TRUE, line=-3)

### END function xplot
}
