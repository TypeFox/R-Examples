# split a 3D object along dimension dim, according to the proportions or
# frequencies specified in vector p

split3d <- function(obj, ...) {
	UseMethod("split3d")
}

split3d.shape3d <- function(obj, p, dim, space=.10, ...) {
	range <-range3d(obj)
	min <- range[1,]
	p <- p/sum(p)                 # assure proportions
	uspace <- space/(length(p)-1) # unit space between objects
	scales <- p * (1-space)
	shifts <- c(0, cumsum(p)[-length(p)])*diff(range[,dim])
	result <- list()
	for (i in seq_along(p)) {
		xscale <- yscale <- zscale <- 1
		xshift <- yshift <- zshift <- 0
		
		if (dim == 1 || tolower(dim)=='x') {
			xscale <- scales[i]
			xshift <- shifts[i] + min[1]*(1-xscale) + (uspace * (i-1))
		} else if (dim == 2|| tolower(dim)=='y') {
			yscale <- scales[i]
			yshift <- shifts[i] + min[2]*(1-yscale) + (uspace * (i-1))
		} else if (dim == 3|| tolower(dim)=='y') {
			zscale <- scales[i]
			zshift <- shifts[i] + min[3]*(1-zscale) + (uspace * (i-1))
		}
		
		result[[i]] <- rgl::translate3d(rgl::scale3d(obj, xscale, yscale, zscale),
				xshift, yshift, zshift)
		
	}
	result
}

# split a list of 3D objects, according to the proportions specified in
# the columns of p.

split3d.list <- function(obj, p, dim, space=.10, ...) {
	nl <- length(obj)
	if (!is.matrix(p) || ncol(p) != nl) stop(gettextf("p must be a matrix with %i columns", nl))
	sl <- list()
	for (i in seq_along(obj)) {
		sl <- c(sl, split3d(obj[[i]], p[,i], dim=dim, space=space))
	}
	sl	
}

#range3d <- function(obj, ...) {
#	UseMethod("range3d")
#}

range3d <- function(obj) {
	if (!"vb" %in% names(obj)) stop("Not a mesh3d or shape3d object")
  x <- with(obj, range(vb[1,]/vb[4,]))
  y <- with(obj, range(vb[2,]/vb[4,]))
  z <- with(obj, range(vb[3,]/vb[4,]))
  result <- cbind(x,y,z)
  rownames(result)<- c('min', 'max')
  result
}

center3d <- function(obj) {
	range <-range3d(obj)
	colMeans(range)
}
