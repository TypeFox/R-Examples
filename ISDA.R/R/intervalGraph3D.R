intervalGraph3D <-
function(xIntervals,yIntervals,zIntervals,yScale){
	
	stopifnot(require("scatterplot3d"))
	## get the vectors of the limits from the Intervals passed as parameter
	xmin = xIntervals$minValue
	xmax = xIntervals$maxValue
	ymin = yIntervals$minValue
	ymax = yIntervals$maxValue
	zmin = zIntervals$minValue
	zmax = zIntervals$maxValue
	## get the min and max limits for x, y and z variables
	minx = min(xmin)
	miny = min(ymin)
	minz = min(zmin)
	maxx = max(xmax)
	maxy = max(ymax)
	maxz = max(zmax)
	
	## drawing plot axes
	z <- seq(minz,maxz,maxz - minz/4 )
	x <- seq(minx,maxx,maxx - minx/4)
	y <- seq(miny,maxy,maxy - miny/4)
	par(xpd = TRUE)
	rr <- scatterplot3d(x,y,z, color = "transparent", box = FALSE, angle = 24,
			xlim = c(minx, maxx), ylim = c(miny, maxy), zlim = c(minz, maxz),scale.y = yScale,grid = FALSE)
	
	cex = 2
	## drawing intervals
	createsRectangles <- function(n,...){
		cubo <- rbind(c(xmin[n],ymin[n],zmax[n]), c(xmin[n],ymin[n],zmin[n]), c(xmax[n],ymin[n],zmin[n]), c(xmax[n],ymax[n],zmin[n]), c(xmax[n],ymax[n],zmax[n]), c(xmin[n],ymax[n],zmax[n]), # < 6 outer
				c(xmax[n],ymin[n],zmax[n]), c(xmin[n],ymax[n],zmin[n]))
		
		## visibile corners + lines:
		rr$points3d(cubo[c(1:6,1,7,3,7,5) ,], cex = cex, type = 'l', lty = 1)
		## hidden corner + lines
		rr$points3d(cubo[c(2,8,4,8,6),     ], cex = cex, type = 'l', lty = 3)
		
	}
	for (i in 1:length(xmin)) {
		
		createsRectangles(i)
	}
}

