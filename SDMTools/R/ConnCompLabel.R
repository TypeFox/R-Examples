#' Connected Components Labelling -- Unique Patch Labelling
#' 
#' \code{ConnCompLabel} is a 1 pass implementation of connected components
#' labelling. Here it is applied to identify disjunt patches within a
#' distribution. \cr \cr The raster matrix can be a raster of class 'asc'
#' (adehabitat package), 'RasterLayer' (raster package) or
#' 'SpatialGridDataFrame' (sp package).
#' 
#' 
#' @param mat is a binary matrix of data with 0 representing background and 1
#' representing environment of interest. NA values are acceptable. The matrix
#' can be a raster of class 'asc' (this & adehabitat package), 'RasterLayer'
#' (raster package) or 'SpatialGridDataFrame' (sp package)
#' @return A matrix of the same dim and class of \code{mat} in which unique
#' components (individual patches) are numbered 1:n with 0 remaining background
#' value.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @seealso \code{\link{PatchStat}}, \code{\link{ClassStat}}
#' @references Chang, F., C.-J. Chen, and C.-J. Lu. 2004. A linear-time
#' component-labeling algorithm using contour tracing technique. Comput. Vis.
#' Image Underst. 93:206-220.
#' @examples
#' 
#' 
#' #define a simple binary matrix
#' tmat = { matrix(c( 0,0,0,1,0,0,1,1,0,1,
#'                    0,0,1,0,1,0,0,0,0,0,
#'                    0,1,NA,1,0,1,0,0,0,1,
#'                    1,0,1,1,1,0,1,0,0,1,
#'                    0,1,0,1,0,1,0,0,0,1,
#'                    0,0,1,0,1,0,0,1,1,0,
#'                    1,0,0,1,0,0,1,0,0,1,
#'                    0,1,0,0,0,1,0,0,0,1,
#'                    0,0,1,1,1,0,0,0,0,1,
#'                    1,1,1,0,0,0,0,0,0,1),nr=10,byrow=TRUE) }
#' 					
#' #do the connected component labelling
#' ccl.mat = ConnCompLabel(tmat)
#' ccl.mat
#' image(t(ccl.mat[10:1,]),col=c('grey',rainbow(length(unique(ccl.mat))-1)))
#' 
#' 
#' @export 
#' @useDynLib SDMTools ccl
ConnCompLabel <- 
function(mat)	{
	#check input for class for returning info
	if (any(class(mat) == 'asc')) { 
		attrib = attributes(mat)
	} else if (any(class(mat) %in% 'RasterLayer')) {
		attrib = mat; mat = asc.from.raster(mat)
	} else if (any(class(mat) == 'SpatialGridDataFrame')) {
		attrib = mat; mat = asc.from.sp(mat)
	} else {
		attrib = attributes(mat)
	}
	#check to ensure matrix
	mat = try(as.matrix(mat))
	if (!is.matrix(mat)) stop('objects must be a matrix')
	#run the connected component labelling
	out = .Call('ccl',mat,PACKAGE='SDMTools')
	#reset the attributes of the input
	if (any(class(attrib) %in% 'RasterLayer')) {
		attrib = setValues(attrib, as.vector(t(t(unclass(out))[dim(out)[2]:1,]))); return(attrib)
	} else if (any(class(attrib) == 'SpatialGridDataFrame')) {
		attrib@data[1] = as.vector(unclass(out)[,dim(out)[2]:1]); return(attrib)
	} else {
		attributes(out) = attrib; return(out)
	}
}
