##' Rotates the x-y coordinates by choosing one datapoint as the origin, choosing another to be fixed on the x axis, and choosing a third to be positive or negative. 
##'
##' Many algorithms exist for projecting m-dimensional datapoints in two-dimensions (e.g., tsne and MDS). 
##' However, often they begin the algorithm by randomly placing datapoints in an arbitrary position. Unfortunately, this makes the axes 
##' meaningless from one iteration of the algorithm to the next, making comparisons across datasets (for example) impossible. One solution 
##' is to fix one datapoint to the origin, while rotating the others about the origin. This algorithm does just that by using polar coordinates.
##'
##' @title Rotate a graph using polar coordinates
##' @aliases rotategraph rotate.graph
##' @param coords a nx2 dimensional matrix corresponding to the cartesian coordinates of the n datapoints.
##' @param scale a nx2 dimensional matrix. This is used when comparing two graphs to one another. This parameter will 
##' scale one graph to the other by ensuring the average of the radi are the same across the two graphs. Defaults to NULL.
##' @param origin an integer indicating which datapoint (i.e., which row) should be fixed as the origin. 
##' @param axis an integer indicating which datapoint (i.e., which row) should be fixed on the x-axis.
##' @param fixedPos an integer indicating which datapoint (i.e., which row) must be positive on y. 
##' @return \item{coords}{the new coordinates obtained after rotation}
##' @return \item{radi}{a vector of the radi for each of the datapoints}
##' @export
##' @seealso \code{\link{compute.theta}}. 
##' @references \url{http://www.mathsisfun.com/polar-cartesian-coordinates.html}
##' @author Dustin Fife
rotateGraph = function(coords, scale=NULL, origin, axis, fixedPos=2){

	#### fix the origin and shift everything else
	new1 = coords - matrix(rep(c(coords[origin,1], coords[origin,2]), times=nrow(coords)), ncol=2, byrow=T)
	
	#### compute polar coordinates
	#### compute distance from "origin" in radius units
	radi = sqrt((coords[,1]-coords[origin,1])^2 + (coords[,2]-coords[origin,2])^2)
	angle = apply(new1, 1, compute.theta)
		
	#### rotate all points the same as B
	angleNew = angle - angle[axis]	
	
	#### scale the radius
	if (!is.null(scale)){
		sc = mean(scale)/mean(radi)
		radi = sc*radi
	}	

	
			#### compute new coordinates
	final = new1
	final[,1] = radi*cos(angleNew)		
	final[,2] = radi*sin(angleNew)		
	final[origin,] = c(0,0)

		#### make sure there's consistency in positive direction
	if (final[fixedPos,2]<0){
		final[,2] = -1*final[,2]
	}		
	
	list(coords = final, radi=radi)
	
}

##' ##' Plots an ellipse
##'
##' @title Plots an Ellipse
##' @param x0 the x coordinate of the center of the circle
##' @param y0 the y coordinate of the center of the circle
##' @param axisX the radius along the x axis
##' @param axisY the radius along the y axis
##' @param color what color should the ellipse be drawn in 
##' @export
##'
##' @author Dustin Fife
ellipse = function(x0, y0, axisX, axisY, color="lightgray"){
	theta <- seq(0, 2 * pi, length=(100))
	x <- x0 + axisX * cos(theta)
	y <- y0 + axisY * sin(theta)
	lines(x, y, type = "l", col=color)
}

##' Computes angle in polar coordinates, accounting for which quadrant the datapoints are in.
##'
##' @title Calculate inverse tangent.
##' @aliases calcT compute.theta
##' @param x a vector of size two that give the cartesian coordinates
##' @return returns the theta, expressed in radians.
##' @author Dustin Fife
##' @export
compute.theta = function(x){
	if ((x[2]>0 & x[1]<0) | (x[2]<0 & x[1]<0)){
		t = atan(x[2]/x[1]) + pi
	} else if (x[2]<0 & x[1]<0){
		t = atan(x[2]/x[1]) + 2*pi
	} else {
		t = atan(x[2]/x[1]) 
	}
	t
}


