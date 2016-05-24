get.box <- function(bathy,x1,x2,y1,y2,width,locator=FALSE,ratio=FALSE,...) {

	if (class(bathy) != "bathy") stop("The matrix provided is not of class bathy")
	if (width<=0) stop("Width must be a positive number")
	if (!locator & (missing(x1)|missing(x2)|missing(y1)|missing(y2))) stop("You need to either use locator=TRUE or specify values for x1, x2, y1 and y2") 

	as.numeric(rownames(bathy)) -> lon
	as.numeric(colnames(bathy)) -> lat

	if (locator) {
		pts <- locator(n=2,type="o",...)
		if (length(pts$x) == 1) stop("Please choose two points from the map")
		x1 <- pts$x[1]
		x2 <- pts$x[2]
		y1 <- pts$y[1]
		y2 <- pts$y[2]
	}

	alpha <- -atan((x2-x1)/(y2-y1))

	p1.x <- x1 + cos(alpha)*width/2
	p2.x <- x2 + cos(alpha)*width/2
	p3.x <- x2 - cos(alpha)*width/2
	p4.x <- x1 - cos(alpha)*width/2

	p1.y <- y1 + sin(alpha)*width/2
	p2.y <- y2 + sin(alpha)*width/2
	p3.y <- y2 - sin(alpha)*width/2
	p4.y <- y1 - sin(alpha)*width/2

	which.min(abs(lon-p1.x)) -> p1x
	which.min(abs(lat-p1.y)) -> p1y
	which.min(abs(lon-p4.x)) -> p4x
	which.min(abs(lat-p4.y)) -> p4y

	which.min(abs(lon-p2.x)) -> p2x
	which.min(abs(lat-p2.y)) -> p2y
	which.min(abs(lon-p3.x)) -> p3x
	which.min(abs(lat-p3.y)) -> p3y

	if (p1x==p4x | p1y==p4y) {
		coord1 <- matrix(as.vector(bathy[p1x:p4x, p1y:p4y]),ncol=length(p1y:p4y),nrow=length(p1x:p4x),dimnames=list(lon[p1x:p4x],lat[p1y:p4y]))
		coord1 <- cbind(as.numeric(dimnames(coord1)[[1]]),as.numeric(dimnames(coord1)[[2]]))
	} else {
		coord1 <- diag.bathy(bathy[p1x:p4x,p1y:p4y],coord=TRUE)
	}

	if (p2x==p3x | p2y==p3y) {
		coord2 <- matrix(as.vector(bathy[p2x:p3x, p2y:p3y]),ncol=length(p2y:p3y),nrow=length(p2x:p3x),dimnames=list(lon[p2x:p3x],lat[p2y:p3y]))
		coord2 <- cbind(as.numeric(dimnames(coord2)[[1]]),as.numeric(dimnames(coord2)[[2]]))
	} else {
		coord2 <- diag.bathy(bathy[p2x:p3x,p2y:p3y],coord=TRUE)
	}

	n1 <- nrow(coord1)
	n2 <- nrow(coord2)
	if (n1<n2) coord2 <- coord2[1:nrow(coord1),]
	if (n1>n2) coord1 <- coord1[1:nrow(coord2),]
	tr <- cbind(coord1,coord2)

	out <- apply(tr,1,function(x) get.transect(x1=x[1],x2=x[3],y1=x[2],y2=x[4],mat=bathy,distance=TRUE))

	di <- round(out[[1]][,3],2)
	prof <- sapply(out,function(x) x[,4])
	
	if (is.list(prof)) {
		nr <- sapply(prof,length)
		prof <- sapply(prof, function(x) x<-x[1:min(nr)])
	}
	
	rownames(prof) <- round(di[1:nrow(prof)])
	colnames(prof) <- round(seq(from=-width/2,to=width/2,len=min(n1,n2)),2)

	deg2km <- function(x1, y1, x2, y2) {

		x1 <- x1*pi/180
		y1 <- y1*pi/180
		x2 <- x2*pi/180
		y2 <- y2*pi/180

		dx <- x2-x1
		dy <- y2-y1
		
		fo <- sin(dy/2)^2 + cos(y1) * cos(y2) * sin(dx/2)^2
		fos <- 2 * asin(min(1,sqrt(fo)))
		
		return(6371 * fos)
	}
	
	d <- max(as.numeric(rownames(prof)))
	pmax <- -min(prof)/1000
	w <- deg2km(p1.x,p1.y,p4.x,p4.y)
	ratios <- round(c(w/d,pmax/d),3)

	if (!locator) # plot(bathy,im=TRUE)
	lines(c(x1,x2),c(y1,y2),type="o",...)
	lines(c(p1.x,p2.x,p3.x,p4.x,p1.x),c(p1.y,p2.y,p3.y,p4.y,p1.y),...)
	
	if (ratio) return(list(depth=prof,ratios=ratios)) else return(prof)
}