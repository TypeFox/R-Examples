get.transect = function(mat, x1, y1, x2, y2, locator=FALSE, distance=FALSE,...){

	as.numeric(rownames(mat)) -> lon
	as.numeric(colnames(mat)) -> lat

## test that the locator input data corresponds to two points
	if (locator) {
		locator(n=2,type="o",...)->coord
		if (length(coord$x) == 1) stop("Please choose two points from the map")
		x1=coord$x[1];x2=coord$x[2];y1=coord$y[1];y2=coord$y[2]
	}

## test that the manual input data is within the bounds of mat
	if ( x1<min(lon) | x1>max(lon) ) stop("x1 not within the longitudinal range of mat")
	if ( x2<min(lon) | x2>max(lon) ) stop("x2 not within the longitudinal range of mat")
	if ( y1<min(lat) | y1>max(lat) ) stop("y1 not within the latitudinal range of mat")
	if ( y2<min(lat) | y2>max(lat) ) stop("y2 not within the latitudinal range of mat")

## reduce mat to the bounds of the input area and switch orientation to get values along the (approximate) diagonal
	which.min(abs(lon-x1)) -> x1b
	which.min(abs(lat-y1)) -> y1b
	which.min(abs(lon-x2)) -> x2b
	which.min(abs(lat-y2)) -> y2b
	
	if (x1b==x2b | y1b==y2b) {
		new.bathy <- matrix(as.vector(mat[x1b:x2b, y1b:y2b]),nrow=length(y1b:y2b),ncol=length(x1b:x2b),dimnames=list(lat[y1b:y2b],lon[x1b:x2b]))
		depth <- as.vector(new.bathy)
	} else {
		new.bathy <- t(mat[x1b:x2b, y1b:y2b])
		depth <- diag.bathy(new.bathy)
	}

## check (and fix if needed) dimentions of new matrix
	as.numeric(colnames(new.bathy)) -> lon
	as.numeric(rownames(new.bathy)) -> lat

	if (length(lon) == 1) lon <- rep(lon,length(lat))
	if (length(lat) == 1) lat <- rep(lat,length(lon))

	if (length(lon)<length(lat)) lon <- seq(lon[1],lon[length(lon)],length.out=length(lat))
	if (length(lon)>length(lat)) lat <- seq(lat[1],lat[length(lat)],length.out=length(lon))

## output format: add distances from x1y1 ?
	if(distance == F) return(data.frame(lon, lat, depth))
	if(distance == T){
		
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
		
		dist.km = NULL
		for(i in 1:length(depth)){
			dist.km = c(dist.km, deg2km(x1=lon[1],y1=lat[1],x2=lon[i],y2=lat[i]))
		}
		out <- data.frame(lon, lat, dist.km, depth)
		return(out)
	}
}