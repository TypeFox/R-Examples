setGeneric("calculate.DVH",
	function (x, dose, ...) {
		standardGeneric("calculate.DVH") 
	}
)

setMethod("calculate.DVH", c("RTdata", "missing"),
	function (x, dose, resolution.xyz=c(0.2,0.2,NA), resolution.dose=0.01, method=NULL, dose.units=NULL) {
		return(calculate.DVH(x$structures, x$dose, resolution.xyz, resolution.dose, method, dose.units))
	}
)

setMethod("calculate.DVH", c("RTdata", "array"),
	function (x, dose, resolution.xyz=c(0.2,0.2,NA), resolution.dose=0.01, method=NULL, dose.units=NULL) {
		return(calculate.DVH(x$structures, dose, resolution.xyz, resolution.dose, method, dose.units))
	}
)

setMethod("calculate.DVH", c("structure3D", "array"),
	function(x, dose, resolution.xyz=c(0.2,0.2,NA), resolution.dose=0.01, method=NULL, dose.units=NULL) {
		method <- match.arg(method, choices=c("ATC", "surface", "axial"))	
		switch(method,
			ATC = return(calc.DVH.ATC(x, dose, resolution.xyz, resolution.dose, dose.units)),
			surface = return(calc.DVH.surface(x, dose, resolution.xyz, resolution.dose, dose.units)),
			axial = return(calc.DVH.axial(x, dose, resolution.xyz, resolution.dose, dose.units)),
			{
				warning("Inappropriate method selection (", method, ")")
				return()
			}
		)
	}
)

setMethod("calculate.DVH", c("structure.list", "array"),
	function(x, dose, resolution.xyz=c(0.2,0.2,NA), resolution.dose=0.01, method=NULL, dose.units=NULL) {
		method <- match.arg(method, choices=c("ATC", "surface", "axial"))	
		switch(method,
			ATC = {
				return(as(lapply(x, function(structure) {
					return(calc.DVH.ATC(structure, dose, resolution.xyz, resolution.dose, dose.units))
				}), "DVH.list"))				
			},
			surface = {
				return(as(lapply(x, function(structure) {
					return(calc.DVH.surface(structure, dose, resolution.xyz, resolution.dose, dose.units))
				}), "DVH.list"))
			},
			axial = {
				return(as(lapply(x, function(structure) {
					return(calc.DVH.axial(structure, dose, resolution.xyz, resolution.dose, dose.units))
				}), "DVH.list"))
			},
			{
				warning("Inappropriate method selection (", method, ")")
				return()
			}
		)
	}
)

setMethod("calculate.DVH", c("ANY", "missing"),
	function (x, dose, ...) {
		warning("Argument 'dose' is missing with no default")
		return()
	}
)

setMethod("calculate.DVH", c("ANY", "array"),
	function (x, dose, ...) {
		warning("Argument 'x' is not an object of class structure3D, structure.list, or RTdata")
		return()
	}
)

setMethod("calculate.DVH", c("ANY", "ANY"),
	function (x, dose, ...) {
		warning("Improper input(s) 'x' and/or 'dose', please refer to RadOnc package documentation for further information")
		return()
	}
)


calc.DVH.ATC <- function(x, dose, resolution.xyz=c(0.2,0.2,NA), resolution.dose=0.01, dose.units=NULL) {
#	reference is Straube and Matthews and Bosche and Prudy 2005 (DVH Analysis: consequences for quality assurance of multi-insitutional clinial trials)
#	placeholder for implementation of multiple DVH calculation algorithms (currrently implemented algorithm makes use of evenly spaced voxel grid and tests every point to see whether or not inside every slice on axial-by-axial basis -- note that this assumes "straight walls" i.e. stepwise surface, no angles/curves between z slices)
#	next step will be to take 3D surface and interpolate dose along z direction as well!!!!
#   also consider different bounding box for each slice?  or use same larger bounding box for overall?  I'll implement both ways and do a timing test to see which one compares better . . . keep current function as-is as baseline check to see performance difference!!!!!!!!!
	if (length(attr(dose, "dose.units")) > 0) {
		dose.units <- attr(dose, "dose.units")	
	}
	dose.units <- match.arg(dose.units, choices=c("cGy", "Gy"))	
	if (length(x$closed.polys) < 1) {
		warning("Structure '", names(x), "' is empty (it contains no pre-defined axial slices)")	
		return()	
	}
	if (dim(x$vertices)[1] <= 2) {
		warning("Structure '", names(x), "' must contain at least three points to calculate a DVH")	
		return()	
	}
	dose.xrange <- range(as.numeric(dimnames(dose)[[1]]), na.rm=TRUE)
	dose.yrange <- range(as.numeric(dimnames(dose)[[2]]), na.rm=TRUE)
	dose.zrange <- range(as.numeric(dimnames(dose)[[3]]), na.rm=TRUE)
	range.struct <- range(x)
	if (any(is.na(range.struct))) {
		warning("Structure '", names(x), "' contains undefined points (coordinate range is undefined)")
		return()
	}
	if ((range.struct[1,1] < dose.xrange[1]) | (range.struct[2,1] > dose.xrange[2]) | 
		(range.struct[1,2] < dose.yrange[1]) | (range.struct[2,2] > dose.yrange[2]) |
		(range.struct[1,3] < dose.zrange[1]) | (range.struct[2,3] > dose.zrange[2])) {
			warning("Structure '", names(x), "' extends beyond calculated dose grid")
			return()
	}
	z.unique <- sort(unique(x$vertices[,3]))
	if (is.na(resolution.xyz[3])) {
		resolution.xyz[3] <- median(abs(diff(z.unique)))
	}
	offset.x <- ((range.struct[2,1]-range.struct[1,1]) %% resolution.xyz[1]) / 2
	if (is.na(offset.x)) {
		warning("Structure '", names(x), "' is not three-dimensional (all points coplanar along x-axis)")
		return()
	}
	xseq <- seq(from=range.struct[1,1]+offset.x, to=range.struct[2,1]-offset.x, by=resolution.xyz[1])
	N.x <- length(xseq)
	offset.y <- ((range.struct[2,2]-range.struct[1,2]) %% resolution.xyz[2]) / 2
	if (is.na(offset.x)) {
		warning("Structure '", names(x), "' is not three-dimensional (all points coplanar along y-axis)")
		return()
	}
	yseq <- seq(from=range.struct[1,2]+offset.y, to=range.struct[2,2]-offset.y, by=resolution.xyz[2])
	N.y <- length(yseq)
	poly.z <- unlist(lapply(x$closed.polys, function(poly) {return(poly[1,3])}))
	N.z <- length(z.unique)
	if (dose.units == "cGy") {
		resolution.dose <- resolution.dose * 100
	}
	doses <- seq(from=min(dose), to=max(dose), by=resolution.dose)
	voxels <- c()
	for (i in z.unique) {
		poly.i <- which(poly.z == i)
		# TEST WHETHER EACH POINT IN GRID IS CONTAINED WITHIN POLYGON(S)
		pts.xyz <- cbind(rep(xseq, each=N.y), rep(yseq, N.x), i)
		results <- rep(0, N.x*N.y)
		lapply(x$closed.polys[poly.i], function(poly) {
			results <<- results+pointInPoly2D(pts.xyz[,1:2], poly[,1:2])
		})
		pts.xyz <- pts.xyz[which(results %%2 != 0),]
		# CALCULATE (INTERPOLATE) DOSES FOR EACH POINT CONTAINED IN STRUCTURE
		dose.xyz <- approx3D(dose, pts.xyz)
		voxels <- c(voxels, dose.xyz)
	}
	dvh <- hist(voxels, breaks=c(doses, max(dose)),plot=FALSE,right=FALSE)$counts*prod(resolution.xyz)/1000
	dvh.volume <- sum(dvh)
	return(new("DVH", type="differential", dose.type="absolute", volume.type="absolute", structure.volume=dvh.volume, doses=doses, volumes=dvh, dose.max=max(voxels,na.rm=TRUE), dose.min=min(voxels,na.rm=TRUE), dose.mean=sum(dvh*doses)/dvh.volume, dose.units=dose.units, structure.name=names(x)))
}


calc.DVH.surface <- function(x, dose, resolution.xyz=c(0.2,0.2,NA), resolution.dose=0.01, dose.units=NULL) {
	if (length(x$triangles) < 1) {
		warning("Structure '", names(x), "' contains no pre-defined triangular mesh")	
		return()	
	}
	if (length(attr(dose, "dose.units")) > 0) {
		dose.units <- attr(dose, "dose.units")	
	}
	dose.units <- match.arg(dose.units, choices=c("cGy", "Gy"))	
	warning("calculate.DVH(..., method='surface') not currently supported")
	return()
}

calc.DVH.axial <- function(x, dose, resolution.xyz=c(0.2,0.2,NA), resolution.dose=0.01, dose.units=NULL) {
#	reference is Straube and Matthews and Bosche and Prudy 2005 (DVH Analysis: consequences for quality assurance of multi-insitutional clinial trials)
#	placeholder for implementation of multiple DVH calculation algorithms (currrently implemented algorithm makes use of evenly spaced voxel grid and tests every point to see whether or not inside every slice on axial-by-axial basis -- note that this assumes "straight walls" i.e. stepwise surface, no angles/curves between z slices)
#	next step will be to take 3D surface and interpolate dose along z direction as well!!!!
#   also consider different bounding box for each slice?  or use same larger bounding box for overall?  I'll implement both ways and do a timing test to see which one compares better . . . keep current function as-is as baseline check to see performance difference!!!!!!!!!
	if (length(attr(dose, "dose.units")) > 0) {
		dose.units <- attr(dose, "dose.units")	
	}
	dose.units <- match.arg(dose.units, choices=c("cGy", "Gy"))	
	if (length(x$closed.polys) < 1) {
		warning("Structure '", names(x), "' is empty (it contains no pre-defined axial slices)")	
		return()	
	}
	dose.xrange <- range(as.numeric(dimnames(dose)[[1]]), na.rm=TRUE)
	dose.yrange <- range(as.numeric(dimnames(dose)[[2]]), na.rm=TRUE)
	dose.zrange <- range(as.numeric(dimnames(dose)[[3]]), na.rm=TRUE)
	range.struct <- range(x)
	if ((range.struct[1,1] < dose.xrange[1]) | (range.struct[2,1] > dose.xrange[2]) | 
		(range.struct[1,2] < dose.yrange[1]) | (range.struct[2,2] > dose.yrange[2]) |
		(range.struct[1,3] < dose.zrange[1]) | (range.struct[2,3] > dose.zrange[2])) {
			warning("Structure '", names(x), "' extends beyond calculated dose grid")
			return()
	}
	z.unique <- sort(unique(x$vertices[,3]))
	if (is.na(resolution.xyz[3])) {
		resolution.xyz[3] <- median(abs(diff(z.unique)))
	}
	offset.x <- ((range.struct[2,1]-range.struct[1,1]) %% resolution.xyz[1]) / 2
	xseq <- seq(from=range.struct[1,1]+offset.x, to=range.struct[2,1]-offset.x, by=resolution.xyz[1])
	N.x <- length(xseq)
	offset.y <- ((range.struct[2,2]-range.struct[1,2]) %% resolution.xyz[2]) / 2
	yseq <- seq(from=range.struct[1,2]+offset.y, to=range.struct[2,2]-offset.y, by=resolution.xyz[2])
	N.y <- length(yseq)
#	offset.z <- ((range.struct[2,3]-range.struct[1,3]) %% resolution.xyz[3]) / 2
#	zseq <- seq(from=range.struct[1,3]+offset.z, to=range.struct[2,3]-offset.z, by=resolution.xyz[3])
#	N.z <- length(zseq)
	poly.z <- unlist(lapply(x$closed.polys, function(poly) {return(poly[1,3])}))
	N.z <- length(z.unique)
	if (dose.units == "cGy") {
		resolution.dose <- resolution.dose * 100
	}
	doses <- seq(from=min(dose), to=max(dose), by=resolution.dose)
	dose.min <- Inf
	dose.max <- -Inf
	dvh.matrix <- matrix(nrow=length(doses),ncol=0)
	for (i in z.unique) {
		poly.i <- which(poly.z == i)
		# TEST WHETHER EACH POINT IN GRID IS CONTAINED WITHIN POLYGON(S)
		pts.xyz <- cbind(rep(xseq, each=N.y), rep(yseq, N.x), i)
		results <- rep(0, N.x*N.y)
		lapply(x$closed.polys[poly.i], function(poly) {
			results <<- results+pointInPoly2D(pts.xyz[,1:2], poly[,1:2])
		})
		pts.xyz <- pts.xyz[which(results %%2 != 0),]
		# CALCULATE (INTERPOLATE) DOSES FOR EACH POINT CONTAINED IN STRUCTURE
		voxels <- approx3D(dose, pts.xyz)
		dose.min <- min(dose.min, min(voxels, na.rm=TRUE), na.rm=TRUE)
		dose.max <- max(dose.max, max(voxels, na.rm=TRUE), na.rm=TRUE)
		dvh.matrix <- cbind(dvh.matrix, hist(voxels, breaks=c(doses, max(dose)),plot=FALSE,right=FALSE)$counts*prod(resolution.xyz)/1000)
	}
	colnames(dvh.matrix) <- z.unique
	dvh.volume <- sum(dvh.matrix)
	return(new("zDVH", type="differential", dose.type="absolute", volume.type="absolute", structure.volume=dvh.volume, doses=doses, volumes=dvh.matrix, dose.max=dose.max, dose.min=dose.min, dose.mean=sum(dvh.matrix*doses)/dvh.volume, dose.units=dose.units, structure.name=names(x)))
}
