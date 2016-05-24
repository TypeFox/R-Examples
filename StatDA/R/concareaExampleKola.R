"concareaExampleKola" <- 
function(x, y, z, zname = deparse(substitute(z)), caname = deparse(substitute(z)), borders="bordersKola", logx = FALSE,
	ifjit = FALSE, ifrev = FALSE, ngrid = 100, ncp = 0, xlim = NULL, xcoord = "Easting", ycoord = 
	"Northing", ifbw = FALSE, x.logfinetick=c(2,5,10),y.logfinetick=c(2,5,10))
{
	# Original wrapper written by Graeme Bonham-Carter, April 2004, to prepare a concentration-
	# area plot as in GEODAS; input data consist of x, y & z, where x & y are coordinates on a 
	# plane, and z is the measured value at point (x,y).  The function uses the interpolation
	# routine in S+ and assumes that area is proportional to the count of grid points.  To be a
	# reasonable model the data points should be 'evenly' spread over the plane.  Interpolated
	# values outside the convex hull of observed data points are set to NA.  The interpolated
	# grid size is computed as (max(x) - min(x))/ngrid, with a default value of 100 for ngrid.
	# The user is prompted for the upper-left and bottom-right corners of a legend panel.
	# If logx = T the data are log-transformed prior to the interpolation step.  If ifjit = T
	# the x and y coordinates are jittered so that no duplicate locations exist, which can
	# cause function interp to fail.  If ifrev= T the empirical concentration-area function is
	# plotted from lowest value to highest.  Triangulation is used if ncp = 0 (default), values
	# of 2 and above result in partial derivatives being used and increased smoothing.  The 
	# plot x-axes are labelled with zname, this can be set to "" and no label is plotted.  The 
	# interpolated 'map' may be titled, often with the text that would be used as the x-axis 
	# label; if no title is required set caname = "".
	#
	# Example: rg.caplot(UTME/1000,UTMN/1000,Cu,zname="",caname="Cu (mg/kg) in O-horizon soil", 
	#  logx=TRUE,ifrev=TRUE,xcoord="Kola Project UTM Easting (km)", 
	#  ycoord="Kola Project UTM Northing (km)")
	#
	# If a "black and white" image is required for monochrome publication set ifbw = T. 
	#

par(mar=c(4,6,4,2))
u <- na.exclude(cbind(x, y, abs(z)))

dx <- (max(u[, 1]) - min(u[, 1]))/ngrid
xo <- seq(from = min(u[, 1]), to = max(u[, 1]), by = dx)
yo <- seq(from = min(u[, 2]), to = max(u[, 2]), by = dx)
zlgnd <- deparse(substitute(z))
if(logx) {
	u[, 3] <- log10(u[, 3])
	zlgnd <- paste("Log10\n", deparse(substitute(z)))
}
#new <- interp.new(u[, 1], u[, 2], u[, 3], xo, yo, duplicate="median",extrap=TRUE)
new <- mba.surf(cbind(u[, 1], u[, 2], u[, 3]), no.X=length(xo), no.Y=length(yo),
       n=1,m=1,extend=TRUE)


if (is.null(borders)){
  #whichdraw <- matrix(as.vector(new$z), nrow=length(xo))
  whichdraw <- matrix(as.vector(new$xyz.est$z), nrow=length(xo))
}
else {
  bord <- get(eval(borders))
  #in.poly=polygrid(new$x,new$y,borders=cbind(bord$x,bord$y),vec.inout=TRUE)
  #whichdraw=matrix(as.vector(new$z)*in.poly$vec.inout, nrow=length(xo))
  in.poly=polygrid(new$xyz.est$x,new$xyz.est$y,borders=cbind(bord$x,bord$y),vec.inout=TRUE)
  whichdraw=matrix(as.vector(new$xyz.est$z)*in.poly$vec.inout, nrow=length(xo))
}

	znew <- whichdraw[in.poly$vec.inout==TRUE]
	if(logx)
		znew <- 10^znew
	xlim <- range(znew)


qpplot.das(u[,3],qdist=qnorm,xlab=zname,xlim=log10(xlim),ylim=qnorm(c(0.0001,0.9999)),
ylab="Cumulative probability [%]", pch=3,cex=0.7, logx=logx,
logfinetick=x.logfinetick,logfinelab=x.logfinetick,line=FALSE,cex.lab=1.2)
title(paste("Original data (n = ",length(u[,3]),")",sep=""))

qpplot.das(log10(znew),qdist=qnorm,xlab=zname,xlim=log10(xlim),ylim=qnorm(c(0.0001,0.9999)),
ylab="Cumulative probability [%]", pch=3,cex=0.7, logx=logx,
logfinetick=x.logfinetick,logfinelab=x.logfinetick,line=FALSE,cex.lab=1.2)
title(paste("Gridded data (n = ",length(znew),")",sep=""))

# plot map

# generate plot with background
im.br=quantile(u[,3],seq(from=0,to=1,by=0.01))
im.col=gray(seq(from=0.1,to=0.9,length=length(im.br)-1))

par(mar=c(1.5,1.5,1.5,1.5))
plot(u[,1],u[,2],frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
#image(new$x,new$y,whichdraw,breaks=im.br,col=im.col, add = TRUE,cex.lab=1.2)
image(new$xyz.est$x,new$xyz.est$y,whichdraw,breaks=im.br,col=im.col, add = TRUE,cex.lab=1.2)
plotbg(map.col=c("gray","gray","gray","gray"),map.lwd=c(1,1,1,1),add.plot=TRUE)

leg.ypos=seq(from=77.7e5,to=78.8e5,length=100)
rect(rep(7.8e5,99),leg.ypos[1:99],rep(8.0e5,99),leg.ypos[2:100],col=im.col,border=FALSE)
rect(7.8e5,leg.ypos[1],8.0e5,leg.ypos[100],border=1)

leg.ypos=seq(from=77.7e5,to=78.8e5,length=5)
text(rep(7.8e5,8),leg.ypos,round(100*c(0,0.25,0.50,0.75,1),2),pos=2,cex=0.8)
im.br=quantile(10^(u[,3]),c(0,0.25,0.50,0.75,1))
text(rep(8.7e5,8),leg.ypos,round(im.br,2),pos=2,cex=0.8)
text(7.4e5,79e5,"Percentile",cex=0.8)
text(8.35e5,79e5,"mg/kg",cex=0.8)
text(7.4e5,79e5,"Percentile",cex=0.8)
text(8.35e5,79e5,"mg/kg",cex=0.8)


# Concentration area plot:

par(mar=c(4,6,4,2))
	conc <- znew[order(znew)]
	cumarea <- seq(1, length(znew))/length(znew) * 100
	if(!ifrev) {
		conc <- rev(conc)
		plot(conc, cumarea, log = "xy", xlab = zname, ylab = 
			"% Cumulative area > values on x-axis", 
            main=paste("Concentration-area plot (n = ",length(conc),")",sep=""),
			xlim = xlim, pch = 3,cex.lab=1.2,xaxt="n",yaxt="n")
	}
	else plot(conc, cumarea, log = "xy", xlab = zname, ylab = 
			"% Cumulative area < values on x-axis", 
            main=paste("Concentration-area plot (n = ",length(conc),")",sep=""),
			xlim = xlim, pch = 3,cex.lab=1.2,xaxt="n",yaxt="n")

axis(1,at=(a<-sort(c((10^(-50:50))%*%t(x.logfinetick)))),labels=a)
axis(2,at=(a<-sort(c((10^(-50:50))%*%t(y.logfinetick)))),labels=a)

# Grid:
abline(v=sort(c((10^(-50:50)%*%t(x.logfinetick)))),lty=3,col=gray(0.5))
abline(h=sort(c((10^(-50:50)%*%t(y.logfinetick)))),lty=3,col=gray(0.5))

	invisible()
}

