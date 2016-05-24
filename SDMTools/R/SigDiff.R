#' Identify Regions of Significant Differences
#' 
#' \code{SigDiff} computes the significance of the pairwise differences
#' relative to the mean and variance of all differences between the two input
#' datasets. This is useful for identifying regions of significant difference
#' between two datasets (e.g., different DEMs (Januchowski et al. 2010) or
#' different species distribution model predictions (Bateman et al 2010)). \cr
#' \cr \code{ImageDiff} is a wrapper to the image.asc command in adehabitat
#' package that uses the result from \code{SigDiff} to create an image mapping
#' the regions of significant differences (positive and negative). \cr \cr
#' \bold{NOTE:} it is assumed the input data are of the same extent and
#' cellsize.
#' 
#' 
#' @param x a vector or matrix of data; the matrix can be of can be a raster of
#' class 'asc' (adehabitat package), 'RasterLayer' (raster package) or
#' 'SpatialGridDataFrame' (sp package)
#' @param y a vector or matrix of data with the same dimensions and class of
#' 'x'
#' @param pattern logical value defining if differences are respective to
#' relative patterning (TRUE) or absolute values (FALSE)
#' @param tasc a matrix of probability values (0 to 1) likely created by
#' \code{SigDiff}; The matrix can be a raster of class 'asc' (adehabitat
#' package), 'RasterLayer' (raster package) or 'SpatialGridDataFrame' (sp
#' package)
#' @param sig.levels the significance levels to define significantly above and
#' below. Default settings represent significance at the 0.05 level
#' @param tcol a set of 3 colors for use in the image to represent
#' significantly lower or greater, and not significant
#' @param ... other graphical parameters defined by image() or plot()
#' @return \code{SigDiff} returns a vector or matrix of the same dimensions and
#' class of the input representing the significance of the pairwise difference
#' relative to the mean and variance of all differences between the two inputs.
#' \cr \cr \code{ImageDiff} returns nothing but creates an image of the areas
#' of significant differences
#' @author Stephanie Januchowski \email{stephierenee@@gmail.com}
#' @references Januchowski, S., Pressey, B., Vanderwal, J. & Edwards, A. (2010)
#' Characterizing errors in topographic models and estimating the financial
#' costs of accuracy. International Journal of Geographical Information
#' Science, In Press. \cr \cr Bateman, B.L., VanDerWal, J., Williams, S.E. &
#' Johnson, C.N. (2010) Inclusion of biotic interactions in species
#' distribution models improves predictions under climate change: the northern
#' bettong Bettongia tropica, its food resources and a competitor. Journal of
#' Biogeography, In Review.
#' @examples
#' #create some simple objects of class 'asc'
#' tasc = as.asc(matrix(1:50,nr=50,nc=50)); print(tasc)
#' #modify the asc objects so that they are slightly different
#' tasc1 = tasc + runif(n = 2500, min = -1, max = 1)
#' tasc2 = tasc + rnorm(n = 2500, mean = 1, sd = 1)
#' 
#' #create graphical representation
#' par(mfrow=c(2,2),mar=c(1,1,4,1))
#' image(tasc1,main='first grid',axes=FALSE)
#' image(tasc2,main='second grid',axes=FALSE)
#' 
#' #get significant difference by spatial patterning
#' out = SigDiff(tasc1,tasc2)
#' ImageDiff(out,main="Pattern Differences",axes=FALSE)
#' 
#' #get significant difference 
#' out = SigDiff(tasc1,tasc2,pattern=FALSE)
#' ImageDiff(out,main="Absolute Differences",axes=FALSE)
#' legend('topleft',legend=c('-ve','ns','+ve'),title='significance',
#'   fill=terrain.colors(3),bg='white')
#' 
#' 
#' @export
SigDiff = function(x,y,pattern=TRUE){
	#check input for class for returning info
	if (class(x) == 'asc') { 
		attrib = attributes(x)
	} else if (any(class(x) %in% 'RasterLayer')) {
		attrib = x; x = asc.from.raster(x); y = asc.from.raster(y)
	} else if (any(class(x) == 'SpatialGridDataFrame')) {
		attrib = x; x = asc.from.sp(x); y = asc.from.sp(y)
	} else {
		attrib = attributes(x)
	}
	
	if(length(which(dim(x)==dim(y)))!=2) stop('asc objects must be of the same extent')#confirm same extents
	pos = which(is.finite(x)) #positions in the data which have a value (are not NA)
	if(pattern) { px = x[pos]/sum(x[pos]) } else { px = x[pos] } #calculate the proportionate value relative to the sum of values across the raster
	if(pattern) { py = y[pos]/sum(y[pos]) } else { py = y[pos] } #calculate the proportionate value relative to the sum of values across the raster
	diff.xy = px-py
	diff.xy=scale(diff.xy) #create z-scores
	t.sig = pnorm(diff.xy) #get the significance values of the z-scores
	out = x; out[pos] = t.sig #create the output ascii grid and write significance values to it
	
	#reset the attributes of the input
	if (any(class(attrib) %in% 'RasterLayer')) {
		attrib = setValues(attrib, as.vector(t(t(unclass(out))[dim(out)[2]:1,]))); return(attrib)
	} else if (any(class(attrib) == 'SpatialGridDataFrame')) {
		attrib@data[1] = as.vector(unclass(out)[,dim(out)[2]:1]); return(attrib)
	} else {
		attributes(out) = attrib; return(out)
	}
}

#' @rdname SigDiff
#' @export
ImageDiff = function(tasc,sig.levels=c(0.025,0.975),tcol=terrain.colors(3),...){
	#check if raster from sp or raster package and convert if necessary
	if (any(class(tasc) %in% 'RasterLayer')) tasc = asc.from.raster(tasc)
	if (any(class(tasc) == 'SpatialGridDataFrame')) tasc = asc.from.sp(tasc)
	tasc[which(is.finite(tasc) & tasc<=sig.levels[1])] = 9
	tasc[which(is.finite(tasc) & tasc>sig.levels[1] & tasc<sig.levels[2])] = 10
	tasc[which(is.finite(tasc) & tasc<=1)] = 11
	image(tasc,col=tcol,zlim=c(9,11),...)
}
