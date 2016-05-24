plot.deldir <- function(x,add=FALSE,wlines=c('both','triang','tess'),
                        wpoints=c('both','real','dummy','none'),
                        number=FALSE,cex=1,nex=1,col=NULL,lty=NULL,
                        pch=NULL,xlim=NULL,ylim=NULL,xlab='x',ylab='y',
                        showrect=FALSE,...)
{
#
# Function plot.deldir to produce a plot of the Delaunay triangulation
# and Dirichlet tesselation of a point set, as produced by the
# function deldir().
#

# Check that x is of the appropriate class.
if(!inherits(x, "deldir")) 
        stop("Argument \"x\" is not of class deldir.\n")

wlines  <- match.arg(wlines)
wpoints <- match.arg(wpoints)

col <- if(is.null(col)) c(1,1,1,1,1) else rep(col,length.out=5)
lty <- if(is.null(lty)) 1:2 else rep(lty,length.out=2)
pch <- if(is.null(pch)) 1:2 else rep(pch,length.out=2)

plot.del <- switch(wlines,both=TRUE,triang=TRUE,tess=FALSE)
plot.dir <- switch(wlines,both=TRUE,triang=FALSE,tess=TRUE)
plot.rl  <- switch(wpoints,both=TRUE,real=TRUE,dummy=FALSE,none=FALSE)
plot.dum <- switch(wpoints,both=TRUE,real=FALSE,dummy=TRUE,none=FALSE)

delsgs <- x$delsgs
dirsgs <- x$dirsgs
n      <- x$n.data
rw     <- x$rw

if(plot.del) {
	x1<-delsgs[,1]
	y1<-delsgs[,2]
	x2<-delsgs[,3]
	y2<-delsgs[,4]
}

if(plot.dir) {
	u1<-dirsgs[,1]
	v1<-dirsgs[,2]
	u2<-dirsgs[,3]
	v2<-dirsgs[,4]
}

X<-x$summary[,"x"]
Y<-x$summary[,"y"]

if(!add) {
	pty.save <- par()$pty
	on.exit(par(pty=pty.save))
	par(pty='s')
	if(is.null(xlim)) xlim <- rw[1:2]
	if(is.null(ylim)) ylim <- rw[3:4]
	plot(0,0,type='n',xlim=xlim,ylim=ylim,
     		xlab=xlab,ylab=ylab,axes=FALSE,...)
	axis(side=1)
	axis(side=2)
}

if(plot.del) segments(x1,y1,x2,y2,col=col[1],lty=lty[1],...)
if(plot.dir) segments(u1,v1,u2,v2,col=col[2],lty=lty[2],...)
if(plot.rl) {
	x.real <- X[1:n]
	y.real <- Y[1:n]
	points(x.real,y.real,pch=pch[1],col=col[3],cex=cex,...)
}
if(plot.dum) {
	x.dumm <- X[-(1:n)]
	y.dumm <- Y[-(1:n)]
	points(x.dumm,y.dumm,pch=pch[2],col=col[4],cex=cex,...)
}
if(number) {
	xoff <- 0.02*diff(range(X))
	yoff <- 0.02*diff(range(Y))
	text(X+xoff,Y+yoff,1:length(X),cex=nex,col=col[5],...)
}
if(showrect) do.call(rect,as.list(x$rw)[c(1,3,2,4)])
invisible()
}
