#'Project plant coordinates onto a viewing plane
#'
#'@description Transform all leaf edge coordinates onto a viewplane coordinate system, with the 'z' axis pointing toward the viewer.
#'
#'
#'@aliases projectplant print.projectedplant3d plot.projectedplant3d
#'@param plant Object of class 'plant3d' (see \code{\link{constructplant}}.
#'@param azimuth,altitude Azimuth and altitude from which plant is viewed.
#'@param x Object of class 'projectedplant3d'.
#'@param silhouette Add a 2D convex hull (see \code{\link{Silhouette}}).
#'@param xlim,ylim Limits for the X and Y axes.
#'@param leaffill If TRUE, fills leaves with green stuff.
#'@param leafcol If leaffill=TRUE, the color of the leaves.
#'@param zerocenter Whether to shift the plant to X=0 and Y=0.
#'@param xlab,ylab Labels for X and Y axes
#'@param \dots Further parameters passed to plot() (but ignored for print).
#'@return Returns an object of class 'projectedplant3d', with components:
#'
#'leaves - A list of matrices with the coordinates of the leaf edges.
#'Each matrix has columns VX,VY and VZ, which are the viewplane coordinates
#'(with the Z axis pointing toward the viewer).
#'
#'viewbound - List of min and max values of coordinates (minx, maxx,
#'miny, etc.).
#'
#'viewangle - Azimuth and altitude used for projection.
#'@author Remko Duursma
#'@seealso \code{\link{STARbar}}
#'@keywords misc
#'@examples
#'
#'
#'# View a plant from above.
#'# The 2D convex hull is also plotted (the 'silhouette'), 
#'# and the area of the hull is printed on the graph.
#'topview <- projectplant(pilularis, altitude=90, azimuth=180)
#'plot(topview, leaffil=TRUE, silhouette=TRUE)
#'
#'@rdname projectplant
#'@export projectplant
projectplant <- function(plant, azimuth, altitude){

	if(class(plant) != "plant3d")stop("Need object of class 'plant3d'\n")

	# Function to get min and max x,y,z values of the transformed plant.
	# This is to place the grid correctly.
	getviewbound <- function(leaflist){
		maxx <- max(sapply(leaflist, function(x)max(x$XYZ[,1])))
		minx <- min(sapply(leaflist, function(x)min(x$XYZ[,1])))
		maxy <- max(sapply(leaflist, function(x)max(x$XYZ[,2])))
		miny <- min(sapply(leaflist, function(x)min(x$XYZ[,2])))
		maxz <- max(sapply(leaflist, function(x)max(x$XYZ[,3])))
		minz <- min(sapply(leaflist, function(x)min(x$XYZ[,3])))
	return(list(minx=minx,maxx=maxx,miny=miny,maxy=maxy,minz=minz,maxz=maxz))
	}

	oldleaves <- plant$leaves
	N <- length(oldleaves)
	p <- makeviewplane(azimuth,altitude)

	# nleafpoints <- nrow(oldleaves[[1]]$XYZ)
	nleafpoints <- sapply(plant$leaves, function(x)nrow(x$XYZ))
	
	xyzs <- do.call("rbind",lapply(oldleaves, function(x)x$XYZ))

	# Get coordinates in viewing plane, using cross-products between current
	# coordinates and the vectors that represent the new system.
	
	#- vector that represents the beam (i.e. z axis in the viewing plane).
	nrm <- p$z 
	
	#- Convert coordinates for all leaf edges.
	vx <- as.vector(xyzs %*% p$x)
	vy <- as.vector(xyzs %*% p$y)
	vv <- as.vector(xyzs %*% nrm)
	
	#- Find acos angle between viewing direction (solar ray) and leaf normal.
	acosangles <- c()
	for(i in 1:N){
		leafn <- plant$leaves[[i]]$leafnormal
		acosangles[i] <- acosangle(leafn, nrm)
	}
	
	#- newleaves : list of matrices with view plane coordinates of leaf edges.
	newxyzs <- cbind(vx,vy,vv)
	colnames(newxyzs) <- c("VX","VY","VZ")
	newleaves <- vector("list",N)
	
	fromi <- cumsum(c(1,nleafpoints[-length(nleafpoints)]))
	toi <- cumsum(nleafpoints)
	
	for(i in 1:N){                           
		newleaves[[i]]$XYZ <- newxyzs[fromi[i]:toi[i],]
		#newxyzs[(1+(i-1)*nleafpoints[i]):(i*nleafpoints[i]) ,]
		newleaves[[i]]$midribpoints <- plant$leaves[[i]]$midribpoints
	}

	l <- list()
	l$leaves <- newleaves
	l$acosangle <- acosangles
	l$viewbound <- getviewbound(newleaves)
	l$viewangle <- list(azimuth=azimuth, altitude=altitude)
	class(l) <- "projectedplant3d"
	
return(l)
}
