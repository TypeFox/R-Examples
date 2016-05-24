
		
gridtrace <- function(obj, npside=100, returnall=TRUE){

	# Get xy bounds of the leaf edges, to place grid
	bounds <- obj$viewbound

	# Place grid
	xwidth <- bounds$maxx - bounds$minx
	ywidth <- bounds$maxy - bounds$miny
	
	# In rare cases, gridpointwidth too close to numerical bounds.
	if(xwidth < 1e-03){
		bounds$maxx <- bounds$maxx + 0.25*ywidth
		bounds$minx <- bounds$minx - 0.25*ywidth
		xwidth <- bounds$maxx - bounds$minx
	}
	if(ywidth < 1e-03){
		bounds$maxy <- bounds$maxy + 0.5*xwidth
		bounds$miny <- bounds$miny - 0.5*xwidth
		ywidth <- bounds$maxy - bounds$miny
	}
	gridpointwidth <- min(xwidth,ywidth) / npside
	
	Px <- sapply(obj$leaves, function(x)x$XYZ[,1])
	Py <- sapply(obj$leaves, function(x)x$XYZ[,2])
	Pz <- sapply(obj$leaves, function(x)x$XYZ[,3])
		
	# Multiple leaf types: kluge to fix uneven number of points per leaf.
	npoints <- sapply(obj$leaves, function(x)nrow(x$XYZ))
	if(length(unique(npoints)) > 1){	
		maxp <- max(npoints)
		shorter <- which(npoints < maxp)
		for(j in shorter){
			n <- npoints[j]
			Px[[j]] <- c(Px[[j]], rep(Px[[j]][n], maxp-n))
			Py[[j]] <- c(Py[[j]], rep(Py[[j]][n], maxp-n))
			Pz[[j]] <- c(Pz[[j]], rep(Pz[[j]][n], maxp-n))
		}
		Px <- do.call(cbind, Px)
		Py <- do.call(cbind, Py)
		Pz <- do.call(cbind, Pz)
	}
		
	
	griddfr <- expand.grid(X=seq(bounds$minx, bounds$maxx, by=gridpointwidth),
							Y=seq(bounds$miny, bounds$maxy, by=gridpointwidth) )
	
	N1 <- nrow(Px)
	N2 <- ncol(Py)
	NS <- nrow(griddfr)

	# Setup for Fortran interface
	xtest <- matrix(griddfr[,1],ncol=1)
	ytest <- matrix(griddfr[,2],ncol=1)
	totarr <- matrix(rep(-999,NS),ncol=1)
	intleaves <- matrix(rep(-999,NS*N2),ncol=N2)
	zinter <- matrix(rep(-999,NS*N2),ncol=N2)
	totonly <- ifelse(returnall, 0L, 1L)
	storage.mode(totarr) <- "integer"
	storage.mode(intleaves) <- "integer"
	storage.mode(xtest) <- storage.mode(ytest) <- "double"

	f <- .Fortran("RAYTRACE", 
			xtest,
			ytest,
			Px,
			Py,
			Pz,
			as.integer(N1),  # number of points per leaf
			as.integer(N2),  # number of leaves
			as.integer(NS),  # number of test points
			totonly,
			totarr,
			intleaves,
			zinter,
			PACKAGE="YplantQMC"
			)
    # !!! Make sure to update when needed !!!
	r <- as.vector(f[[10]])
	intleaves <- f[[11]]
	zinter <- f[[12]]
	
	# Number of rays that intersect with at least one leaf on the plant.
	nintersect <- length(r[r >= 1])
	
	props <- list(nintersect=nintersect, pointarea=gridpointwidth^2, gridpointwidth=gridpointwidth, npside=npside)
	
	if(!returnall)
		return(props)
	else {
		griddfr$p <- r
		l <- list()
		l$rays <- griddfr
		l$intleaves <- intleaves
		l$zinter <- zinter
		l <- c(l, props)
		class(l) <- "tracedplant"
		return(l)
	}
}


# OLD slow R-native version.
# # Returns a dataframe of x,y coordinates that intersect with the current viewpoint.
# # Should be much less slow for large plants (i.e., computing time scales with much less than
# # nleaves).
# raytracetotal <- function(obj, npside=100){

	# if(class(obj) != "projectedplant3d")stop("Need projected plant.\n")

	# # Get xy bounds of the leaf edges, to place grid
	# bounds <- obj$viewbound

	# # Place grid
	# xwidth <- bounds$maxx - bounds$minx
	# ywidth <- bounds$maxy - bounds$miny
	# gridpointwidth <- min(xwidth,ywidth) / npside
	
	# gridcoor <- expand.grid(x=seq(bounds$minx, bounds$maxx, by=gridpointwidth),
							# y=seq(bounds$miny, bounds$maxy, by=gridpointwidth) )

	# nleaves <- length(obj$leaves)
							
							
							
	# # For all grid points, find if they intersect with leaves.
	# # Function 'point.in.polygon' is from the 'sp' package.
	# thesedontintersect <- 1:nrow(gridcoor)
	# for(i in 1:nleaves){
		# val <- point.in.polygon(gridcoor$x[thesedontintersect], gridcoor$y[thesedontintersect], 
			# obj$leaves[[i]]$XYZ[,1], obj$leaves[[i]]$XYZ[,2])
		
		# # Cull points for testing that have been tested already:
		# thesedontintersect <- thesedontintersect[val==0]
	# }
	# S <- gridcoor[setdiff(1:nrow(gridcoor), thesedontintersect),]
	
# return(list(nintersect=nrow(S), pointarea=gridpointwidth^2))
# }


