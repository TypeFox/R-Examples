# Author: Robert J. Hijmans
# Date :  July 2010
# Version 1.0
# Licence GPL v3

# Based on functions in R package 'RgoogleMaps' 
# by Markus Loecher, Sense Networks <markus at sensenetworks.com>

# October 2012
# Updated with contributions by Sebastien Rochette

gmap <- function(x, exp=1, type='terrain', filename='', style=NULL, scale=1, zoom=NULL, size=c(640, 640), rgb=FALSE, lonlat=FALSE, ...) {

	if (! requireNamespace('rgdal') ) { 
		stop('rgdal not available') 
	}
	
	if (! type %in% c('roadmap', 'satellite', 'hybrid', 'terrain')) {
		warning("type should be: roadmap, satellite, hybrid, or terrain: Terrain chosen by default") 
		type <- 'terrain'
	}
	
	mxzoom <- function(latrange, lonrange, size=size) {
        SinPhi = sin(latrange * pi/180)
        normX = lonrange/180
        normY = (0.5 * log(abs((1 + SinPhi)/(1 - SinPhi))))/pi
        MaxZoom.lon <- floor(1 + log2(abs(size[1]/256/diff(normX))))
        MaxZoom.lat <- floor(1 + log2(abs(size[2]/256/diff(normY))))
        return(c(MaxZoom.lat = MaxZoom.lat, MaxZoom.lon = MaxZoom.lon))
    }	

	ll2XY <- function (lat, lon, zoom) {
	# function from in R package 'RgoogleMaps' 
	# by Markus Loecher, Sense Networks <markus at sensenetworks.com>
		SinPhi = sin(lat * pi/180)
		normX = lon/180
		normY = (0.5 * log((1 + SinPhi)/(1 - SinPhi)))/pi
		Y = (2^zoom) * ((1 - normY)/2)
		X = (2^zoom) * ((normX + 1)/2)
		x = 256 * (X - floor(X))
		y = 256 * (Y - floor(Y))
		return(list(Tile = cbind(X = floor(X), Y = floor(Y)), Coords = cbind(x = x,  y = y)))
	}

	xy2ll <- function (MyMap, X, Y) {
	# function from in R package 'RgoogleMaps' 
	# by Markus Loecher, Sense Networks <markus at sensenetworks.com>
		lat.center <- MyMap[[1]]
		lon.center <- MyMap[[2]]
		zoom <- MyMap[[3]]
		mycenter <- ll2XY(lat.center, lon.center, zoom)
		x <- mycenter$Tile[, "X"] + (X + mycenter$Coords[, "x"])/256
		y <- mycenter$Tile[, "Y"] - (Y - mycenter$Coords[, "y"])/256
		ytilde <- 1 - y/2^(zoom - 1)
		yy = (exp(2 * pi * ytilde) - 1)/(exp(2 * pi * ytilde) + 1)
		ShiftLat <- function(yy) {
			n = c(-1, 0, 1)
			lat = 2 * pi * (n) + asin(yy)
			lat <- lat[which(lat <= pi/2 & lat > -pi/2)]
			lat <- 180 * lat/pi
			return(lat)
		}
		lat <- sapply(yy, ShiftLat)
		lon = 180 * (x/2^(zoom - 1) - 1)
		return(cbind(lat = lat, lon = lon))
	}

	tile2r <- function (points, center) {
	# function from in R package 'RgoogleMaps' 
	# by Markus Loecher, Sense Networks <markus at sensenetworks.com>
		X <- 256 * (points$Tile[, "X"] - center$Tile[, "X"]) + (points$Coords[, "x"] - center$Coords[, "x"])
		Y <- -256 * (points$Tile[, "Y"] - center$Tile[, "Y"]) - (points$Coords[, "y"] - center$Coords[, "y"])
		return(list(X = X, Y = Y))
	}

#	gurl <- "http://maps.google.com/staticmap?"
	gurl <- "http://maps.googleapis.com/maps/api/staticmap?"
	
	if (is.character(x)) {
		x <- geocode(x, oneRecord=TRUE)
		if (is.na(x$latitude)) { 
			stop('location not found') 
		}
		x <- extent(as.vector(as.matrix(x[5:8])))

	} else {
		prj <- projection(x, asText=TRUE)
		
		if ( isLonLat(prj) ) {
		
			x <- extent(x)
		} else {
			if ( is.na(prj) ) {
				bb <- extent(x)
				extLL <- (bb@xmin > -366 & bb@xmax < 366 & bb@ymin > -90.1 & bb@ymax < 90.1) 
				if (extLL) {
					x <- bb
				} else {
				# warning('CRS is unknown, and does not look like Lon/Lat, assuming it is Mercator')
					rad <- 6378137	
					p <- t(bbox(x)) 
					p[, 2] <- pi/2 - 2 * atan(exp(-p[, 2]/rad))
					p[, 1] <- p[, 1]/rad
					p <- p / (pi/180)
					x <- extent(p[1,1], p[2,1], p[1,2], p[2,2])
				}
			} else {
				x <- extent( projectExtent(x, "+proj=longlat +datum=WGS84") )
			}
		}
	}
		e <- x * exp
		e@xmin <- max(-180, e@xmin)
		e@xmax <- min(180, e@xmax)
		e@ymax <- min(89, e@ymax)
		e@ymin <- max(-89, e@ymin)
		
		lonR <- c(e@xmin, e@xmax)
		latR <- c(e@ymin, e@ymax)

# 		size <- c(640, 640)
		if(is.null(zoom)){
		  zoom <- min(mxzoom(latR, lonR, size))
		}
		center <- c(mean(latR), mean(lonR))
 	
		ll <- ll2XY(latR[1], lonR[1], zoom)
		ur <- ll2XY(latR[2], lonR[2], zoom)
		cr <- ll2XY(center[1], center[2], zoom)
		ll.Rcoords <- tile2r(ll, cr)
		ur.Rcoords <- tile2r(ur, cr)
		s1 <- 2 * max(c(ceiling(abs(ll.Rcoords$X)), ceiling(abs(ur.Rcoords$X)))) +   1
		s2 <- 2 * max(c(ceiling(abs(ll.Rcoords$Y)), ceiling(abs(ur.Rcoords$Y)))) +   1
		if (s1 > s2) {
			size[2] <- as.integer(size[2] * s2/s1)
		} else {
			size[1] <- as.integer(size[1] * s1/s2)
		}
		s <- paste(size, collapse = "x")

		ctr <- paste(center, collapse = ",")
	
		gurl <- paste(gurl, "center=", ctr, "&zoom=", zoom, "&size=", s, "&maptype=", type, "&format=gif&sensor=false&scale=", scale, sep = "")
		if (!is.null(style)) {
			style <- gsub("\\|", "%7C", style)
			style <- gsub(" ", "", style)
			gurl <- paste(gurl, '&style=', style, sep='')
		}
#		message(gurl, "\n")
	
	filename <- trim(filename)
	if (filename == '') {
		filename <- rasterTmpFile()
	}
	extension(filename) <- 'gif'
	download.file(gurl, filename, mode="wb", quiet=TRUE)
    
	MyMap <- list(lat.center = center[1], lon.center = center[2], zoom = zoom)
	bb <- list(ll = xy2ll(MyMap, X=-size[1]/2 + 0.5, Y=-size[2]/2 - 0.5), ur=xy2ll(MyMap, X=size[1]/2 + 0.5, Y=size[2]/2 - 0.5))

	r <- raster(filename, warn=FALSE)
	ext <- extent(bb$ll[2], bb$ur[2], bb$ll[1], bb$ur[1])
	p <- t(bbox(raster(ext))) *  pi/180
	rad <- 6378137	
    p[,2] <- log(tan(p[, 2]) + (1/cos(p[, 2]))) 
    p <- p * rad
	extent(r) <- extent(as.vector(p))
	projection(r) <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"
	try( hdr(r, format='worldfile', extension='.gfw') )
	extension(filename) <- 'prj'
	rgdal::showWKT(projection(r), file=filename, morphToESRI=TRUE)
	
	if (lonlat) {
		ct <- r@legend@colortable 
		r <- projectRaster(r, crs="+proj=longlat +datum=WGS84", method='ngb')
		r@legend@colortable <- ct
	}
	
    if (rgb) {
		d <- t( col2rgb(r@legend@colortable) )
		d <- data.frame(id=0:255, d)
		r <- subs(r, d, which=2:4)
    }
	
    return(r)
}

#e = extent( -121.9531 , -120.3897 , 35.36 , 36.61956 )
#r = gmap(e)
#plot(r)