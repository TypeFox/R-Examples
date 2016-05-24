dumpts <- function(x,y,dpl,rw) {
#
# Function dumpts to append a sequence of dummy points to the
# data points.
#

ndm   <- 0
xd    <- NULL
yd    <- NULL
xmin  <- rw[1]
xmax  <- rw[2]
ymin  <- rw[3]
ymax  <- rw[4]

# Points on radii of circles emanating from data points:
if(!is.null(dpl$nrad)) {
	nrad  <- dpl$nrad # Number of radii from each data point.
	nper  <- dpl$nper # Number of dummy points per radius.
	fctr  <- dpl$fctr # Length of each radius = fctr * mean
                             # interpoint distance.
	lrad  <- fctr*mnnd(x,y)/nper
	theta <- 2*pi*(1:nrad)/nrad
	cs    <- cos(theta)
	sn    <- sin(theta)
	xt    <- c(lrad*(1:nper)%o%cs)
	yt    <- c(lrad*(1:nper)%o%sn)
	xd    <- c(outer(x,xt,'+'))
	yd    <- c(outer(y,yt,'+'))
}

# Ad hoc points passed over as part of dpl:
if(!is.null(dpl$x)) {
	xd <- c(xd,dpl$x)
	yd <- c(yd,dpl$y)
}

# Delete dummy points outside the rectangular window.
ndm  <- length(xd)
if(ndm >0) {
	drop <- (1:ndm)[xd<xmin|xd>xmax|yd<ymin|yd>ymax]
	if(length(drop)>0) {
		xd  <- xd[-drop]
		yd  <- yd[-drop]
	}
}

# Rectangular grid:
ndx <- dpl$ndx
okx <- !is.null(ndx) && ndx > 0
ndy <- dpl$ndy
oky <- !is.null(ndy) && ndy > 0
if(okx & oky) {
	xt  <- if(ndx>1) seq(xmin,xmax,length=ndx) else 0.5*(xmin+xmax)
	yt  <- if(ndy>1) seq(ymin,ymax,length=ndy) else 0.5*(ymin+ymax)
	xy <- expand.grid(x=xt,y=yt)
	xd  <- c(xd,xy$x)
	yd  <- c(yd,xy$y)
}

ndm <- length(xd)
list(x=c(x,xd),y=c(y,yd),ndm=ndm)
}
