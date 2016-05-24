"hcube"<-
function(x, scale, translation)
{
#   DATE WRITTEN:  24 April 1995          LAST REVISED:  1 May 1995
#   AUTHOR:  Scott D. Chasalow
#
#   DESCRIPTION:
#         Generate all points on a hypercuboid lattice. 
#         Argument x is an integer vector giving the extent of each dimension; 
#         the number of dimensions is length(x).  
#         Argument scale is a vector of real numbers giving an amount by which 
#         to multiply the points in each dimension;  it will be replicated as 
#         necessary to have the same length as x.
#         Argument translate is a vector of real numbers giving an amount to 
#         translate (from the "origin", rep(1,length(x))) the points in each 
#         dimension;  it will be replicated as necessary to have the same 
#         length as x.  To use rep(0,length(x)) as the origin,  use 
#         translation = -1.  Scaling,  if any,  is done BEFORE translation.
#
#   VALUE:
#         A prod(x) by length(x) numeric matrix;  element (i,j) gives the 
#         location of point i in the jth dimension.  The first column 
#         (dimension) varies most rapidly.
#
#   SEE ALSO:
#         fac.design,  expand.grid
#
	ncols <- length(x)
	nrows <- prod(x)
	cp <- c(1, cumprod(x)[ - ncols])
	out <- lapply(as.list(1:length(x)), function(i, a, b, nr)
	rep(rep(1:a[i], rep(b[i], a[i])), length = nr), a = x, b = cp, nr = 
		nrows)
	out <- array(unlist(out), c(nrows, ncols))
	if(!missing(scale)) {
		scale <- rep(scale, length = ncols)
		out <- sweep(out, 2, scale, FUN = "*")
	}
	if(!missing(translation)) {
		translation <- rep(translation, length = ncols)
		out <- sweep(out, 2, translation, FUN = "+")
	}
	out
}

