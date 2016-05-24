#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Douglas-Peucker Polyline Simplification
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copied from shapefiles/R/shapefiles.R 091119 RSB
#Ben Stabler, bstabler@ptvamerica.com, March 2005
#Douglas, D. and Peucker, T. (1973). "Algorithms for 
#the reduction of the number of points required to 
#represent a digitized line or its caricature." 
#The Canadian Cartographer 10(2). 112-122.

#Currently uses the line, not the line segment to 
#determine the distance of the points from the line.
#See http://www.lgc.com/resources/Doug_Peucker.pdf
#for more information.  This can result in the 
#omission of extreme "outlier-like" points.
####################################################

dp_func <- function(points, tolerance) {

 #Convert to lowercase
 names(points) <- tolower(names(points))

 #Calculate distance between two points
 distance <- function(x1,x2,y1,y2) {
 	sqrt((x2-x1)^2 + (y2-y1)^2)
 }
 
 #Calculate equation of a line from two points
 equationOfLine <- function(x1,x2,y1,y2) {
 	slope <- (y2-y1)/(x2-x1)
 	b <- y1 - slope * x1
 	c(slope,b)
 }
 
 #Calculate y-intercept from a point and a slope
 calcB <- function(slope, px, py) {
 	py - (slope * px) 
 }
 
 #Calculate intercept of two lines from their slope and y-intercept
 intercept <- function(s1,b1,s2,b2) {
 	x1 <- (b2-b1)/(s1-s2)
 	y1 <- s1 * x1 + b1
 	c(x1,y1)
 }

 #Setup vector to mark points to keep
 keep <- rep(FALSE, length(points$x))
 keep[1] <- TRUE
 keep[length(keep)] <- TRUE

 #Function definition to simplify points
 simplify <- function(start, end, tol=tolerance) {

  #Calculate intermediate point distances 
  if (length(points$x[start:end]) > 2L) {

 	#Avoid Inf slope
 	if( points$x[start] ==  points$x[end] )  
            points$x[start] <- points$x[start] - 0.0000001
 	if( points$y[start] ==  points$y[end])  
            points$y[start] <- points$y[start] - 0.0000001
 
	#Calculate equation of line of middle points
 	line <- equationOfLine( points$x[start], points$x[end],
            points$y[start], points$y[end])
	
	#Calculate y-intercepts
 	b <- mapply(function(x,y) calcB(-1/line[1], x, y),
            points$x[(start+1):(end-1)], points$y[(start+1):(end-1)])
	
	#Calculate intercepts with with start-end line
 	ints <- sapply(b, function(x)
            intercept(line[1], line[2], -1/line[1], x), simplify=FALSE)

	#Calculate distances of points from line
 	distances <- mapply(function(x,y,z) distance(x[[1]], y, x[[2]], z),
            ints, points$x[(start+1):(end-1)], points$y[(start+1):(end-1)])

	#If any point greater than tolerance split at max distance point 
        #and apply to each side
 	if (any(distances >= tol)) {
		keep[which.max(distances)+start] <<- TRUE
		#print(which.max(distances)+start)
		
		simplify(start, which.max(distances)+start)
		simplify(which.max(distances)+start, end)
	}

  }

 }

 #Start simplification 
 simplify(1, length(points$x))

 #Return simplified points
 list(x=points$x[keep], y=points$y[keep])
}


thinnedSpatialPoly <- function(SP, tolerance, minarea=0, topologyPreserve=FALSE, avoidGEOS=FALSE) {
    stopifnot(inherits(SP, "SpatialPolygons"))
    rgeosI <- rgeosStatus()
    if (rgeosI && !avoidGEOS) {
        # require(rgeos)
    	if (!requireNamespace("rgeos", quietly = TRUE))
			stop("package rgeos required for tinnedSpatialPoly")
        res <- rgeos::gSimplify(spgeom=SP, tol=tolerance, topologyPreserve=topologyPreserve)
    } else {
      pls <- slot(SP, "polygons")
      pls_dp <- vector(mode="list", length=length(pls))
      for (i in 1:length(pls)) {
        Pls <- slot(pls[[i]], "Polygons")
        Pls_dp <- vector(mode="list", length=length(Pls))
        for (j in 1:length(Pls)) {
            crds <- slot(Pls[[j]], "coords")
            crds_s <- dp_func(list(x=crds[,1], y=crds[,2]), tolerance=tolerance)
            crds_s <- do.call("cbind", crds_s)
            if(!identical(crds_s[1,], crds_s[nrow(crds_s),]))
            crds_s <- rbind(crds_s, crds_s[1,])
            Pls_dp[[j]] <- Polygon(crds_s)
        }
        areas <- sapply(Pls_dp, slot, "area")
        Keep <- areas > minarea  
        if (all(!Keep)) Keep[which.max(areas)] <- TRUE
        Pls_dp <- Pls_dp[Keep]
        pls_dp[[i]] <- Polygons(Pls_dp, ID=slot(pls[[i]], "ID"))
      }
      res <- SpatialPolygons(pls_dp, proj4string=CRS(proj4string(SP)))
    }
    if (is(SP, "SpatialPolygonsDataFrame"))
        res <- SpatialPolygonsDataFrame(res, data=slot(SP, "data"))
    res
}

