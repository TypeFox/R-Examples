#' I Similarity Statistic for Quantifying Niche Overlap
#' 
#' \code{Istat} computes the I similarity statistic of Warren et al. 2008. It
#' is a method for defining niche overlap from predictions of species'
#' distributions. \cr \cr \bold{NOTE:} it is assumed the input data are of the
#' same extent and cellsize, and all values are positive.
#' 
#' The I similarity statistic sums the pair-wise differences between two
#' predictions to create a single value representing the similarity of the two
#' distributions. The I similarity statistic ranges from a value of 0, where
#' two distributions have no overlap, to 1 where two distributions are
#' identical (Warren et al., 2008).
#'
#' NOTE: updated to correct equation but not to worry about old... see explanation at
#' \url{http://enmtools.blogspot.com.au/2010_09_01_archive.html}.
#' 
#' @param x a vector or matrix of data; the matrix can be a raster of class
#' 'asc' (adehabitat package), 'RasterLayer' (raster package) or
#' 'SpatialGridDataFrame' (sp package)
#' @param y a vector or matrix of data with the same dimensions and class of
#' 'x'
#' @param old a boolean identifying if "old" equation is to be used (see description).
#' This was kept for legacy issues.
#'
#' @return A single value that is the I similarity statistic
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @references Warren, D. L., R. E. Glor, M. Turelli, and D. Funk. 2008.
#' Environmental Niche Equivalency versus Conservatism: Quantitative Approaches
#' to Niche Evolution. Evolution 62:2868-2883.
#' @examples
#' 
#' 
#' #create some simple objects of class 'asc'
#' tasc = as.asc(matrix(1:50,nr=50,nc=50)); print(tasc)
#' #modify the asc objects so that they are slightly different
#' tasc1 = tasc + runif(n = 2500, min = -1, max = 1)
#' tasc2 = tasc + rnorm(n = 2500, mean = 1, sd = 1)
#' 
#' #ensure all data is positive
#' tasc1 = abs(tasc1)
#' tasc2 = abs(tasc2)
#' 
#' #calculate the I similarity statistic
#' I = Istat(tasc1,tasc2)
#' print(I) #high niche overlap
#' 
#' #using a more variable map
#' tasc2 = tasc + rnorm(n = 2500, mean = 25, sd = 15);tasc2 = abs(tasc2)
#' I = Istat(tasc1,tasc2)
#' print(I) #lower niche overlap
#' 
#' 
#' @export 
Istat = function(x,y,old=FALSE){
	#check if raster from sp or raster package and convert if necessary
	if (any(class(x) %in% 'RasterLayer')) { x = asc.from.raster(x); y = asc.from.raster(y) }
	if (any(class(x) == 'SpatialGridDataFrame')) { x = asc.from.sp(x); y = asc.from.sp(y) }
	if(length(which(dim(x)==dim(y)))!=2) stop('matrix / raster objects must be of the same extent')#confirm same extents
	if (min(c(x,y),na.rm=T)<0) stop('all values must be positive')#confirm all values are positive
	pos = which(is.finite(x) & is.finite(y)) #positions in the data which have a value (are not NA)
	px = x[pos]/sum(x[pos]) #calculate the proportionate value relative to the sum of values across the raster
	py = y[pos]/sum(y[pos]) #calculate the proportionate value relative to the sum of values across the raster
	H = sqrt(sum((sqrt(px)-sqrt(py))^2)) #define the Hellenger distance
	if (old) {
		return (1 - 0.5 * H) #return the old I-Statistic
	} else {
		return (1 - (H^2)/2) #return the corrected I-Statistic
	}
}
