approx3D <- function(data, x, y=NULL, z=NULL, method=c("trilinear"), extrapolate=FALSE) {
	method <- match.arg(method)
	if ((length(x) == 3) & is.null(y) & is.null(z)) {
		z <- as.numeric(x[3])
		y <- as.numeric(x[2])
		x <- as.numeric(x[1])
	}
	if (is.matrix(x) & is.null(y) & is.null(z)) {
		if (dim(x)[2] == 3) {
			z <- as.numeric(x[,3])
			y <- as.numeric(x[,2])
			x <- as.numeric(x[,1])
		}
	}
	if (length(x) != length(y)) {
		stop("'x' and 'y' lengths differ")
	}
	if (length(x) != length(z)) {
		stop("'x' and 'z' lengths differ")
	}
	if (length(x) < 1) {
		return()
	}
	x.coord <- as.numeric(dimnames(data)[[1]])
	y.coord <- as.numeric(dimnames(data)[[2]])
	z.coord <- as.numeric(dimnames(data)[[3]])
	## FUNCTION ASSUMES THAT GRID SPACING IS EVEN IN THE X, Y, AND Z DIRECTIONS
	x.diff <- abs(x.coord[2]-x.coord[1])
	y.diff <- abs(y.coord[2]-y.coord[1])
	z.diff <- abs(z.coord[2]-z.coord[1])
	if (!extrapolate) {
		x[x < min(x.coord)] <- NA
		x[x > max(x.coord)] <- NA
		y[y < min(y.coord)] <- NA
		y[y > max(y.coord)] <- NA
		z[z < min(z.coord)] <- NA
		z[z > max(z.coord)] <- NA
	}
	else {
		x <- pmax(x, min(x.coord))
		x <- pmin(x, max(x.coord))
		y <- pmax(y, min(y.coord))
		y <- pmin(y, max(y.coord))
		z <- pmax(z, min(z.coord))
		z <- pmin(z, max(z.coord))
	}
	x.unique <- unique(x)
	x.unique.data <- lapply(x.unique,
		function(x) {
			if (x %in% x.coord) {
				return(list(coords=rep(which(x.coord == x), 2), dist=0))
			}
			else {
				coords <- order(abs(x.coord-x))[1:2]
				return(list(coords=coords, dist=abs((x-x.coord[coords[1]])/x.diff)))
			}
		}
	)
	y.unique <- unique(y)
	y.unique.data <- lapply(y.unique,
		function(y) {
			if (y %in% y.coord) {
				return(list(coords=rep(which(y.coord == y), 2), dist=0))
			}
			else {
				coords <- order(abs(y.coord-y))[1:2]
				return(list(coords=coords, dist=abs((y-y.coord[coords[1]])/y.diff)))
			}
		}
	)
	names(y.unique.data) <- y.unique
	z.unique <- unique(z)
	z.unique.data <- lapply(z.unique,
		function(z) {
			if (z %in% z.coord) {
				return(list(coords=rep(which(z.coord == z), 2), dist=0))
			}
			else {
				coords <- order(abs(z.coord-z))[1:2]
				return(list(coords=coords, dist=abs((z-z.coord[coords[1]])/z.diff)))
			}
		}
	)
	# TRILINEAR INTERPOLATION METHOD
	if (method == "trilinear") {
		data <- mapply(function(x, y, z) {
				x.i <- x.unique.data[[which(x.unique == x)]]
				y.i <- y.unique.data[[which(y.unique == y)]]
				z.i <- z.unique.data[[which(z.unique == z)]]
				data.xyz <- data[x.i$coords, y.i$coords, z.i$coords]
				data.xyz <- data.xyz[2,,]*x.i$dist + data.xyz[1,,]*(1-x.i$dist)
				data.xyz <- data.xyz[2,]*y.i$dist + data.xyz[1,]*(1-y.i$dist)
				return(as.numeric(data.xyz[2]*z.i$dist + data.xyz[1]*(1-z.i$dist)))
			},
			x, y, z
		)
		return(unlist(data))
	}
	else {
		warning("Only trilinear interpolation method currently supported")
		return(NA)
	}
} 