#' Centre of Gravity or Mass calculations for spatial data
#' 
#' \code{COGravity} calculates the Centre of Gravity (or also known as Centre
#' of Mass) for point or raster spatial data.\cr \cr \bold{Note:} NA data is
#' automatically ommitted from analysis.
#' 
#' For raster-based data, if \code{wt} is missing, the values of the ascii are
#' assumed to be the weights; otherwise, the values are assumed to be the
#' \code{z} values.
#' 
#' @param x a vector of e.g., longitudes or eastings, or a raster of class
#' 'asc', 'RasterLayer' or 'SpatialGridDataFrame'.
#' @param y a vector of e.g., latitudes or northings.
#' @param z a vector of e.g., elevations.
#' @param wt a vector or raster of class 'asc', 'RasterLayer' or
#' 'SpatialGridDataFrame' representing weights for data.
#' @return Returns a named vector of data representing the Centre of Gravity in
#' x, y & z dimensions (depending on data supplied).
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @examples
#' 
#' 
#' #create some points
#' x = seq(154,110,length=25)
#' y = seq(-10,-54,length=25)
#' z = seq(100,200,length=25)
#' wt = runif(25) #random weights
#' #calculate the Centre of Gravity for these points
#' COGravity(x,y,z,wt)
#' 
#' #create a simple objects of class 'asc'
#' x = as.asc(matrix(1:50,nr=50,nc=50))
#' wt = as.asc(matrix(runif(50),nr=50,nc=50))
#' 
#' #calculate COG with weighting defined in x 
#' COGravity(x)
#' #calculate COG with weighting defined in wt (values in x are assumed elevation (z)) 
#' COGravity(x,wt=wt)
#' 
#' 
#' @export COGravity
COGravity <- function(x,y=NULL,z=NULL,wt=NULL) {
	#check if raster from sp or raster package and convert if necessary
	if (any(class(x) %in% 'RasterLayer')) x = asc.from.raster(x)
	if (any(class(x) == 'SpatialGridDataFrame')) x = asc.from.sp(x)
	#check if x is vector or matrix
	if (is.vector(x)) { #if the data is a vector...do calculations
		if (is.null(wt)) {	#if no weighting supplied, calculate means & standard deviations
			out = c(COGx=mean(x,na.rm=TRUE),COGx.sd=sd(x,na.rm=TRUE))
			if (!is.null(y)) out = c(out,COGy=mean(y,na.rm=TRUE),COGy.sd=sd(y,na.rm=TRUE))
			if (!is.null(z)) out = c(out,COGz=mean(z,na.rm=TRUE),COGz.sd=sd(z,na.rm=TRUE))
		} else { #if weighting supplied, calculate weighted means and variances to get COG
			out = c(COGx=wt.mean(x,wt),COGx.sd=wt.sd(x,wt))
			if (!is.null(y)) out = c(out,COGy=wt.mean(y,wt),COGy.sd=wt.sd(y,wt))
			if (!is.null(z)) out = c(out,COGz=wt.mean(z,wt),COGz.sd=wt.sd(z,wt))
		}
	} else if (any(class(x) == 'asc')) { #if x is of class 'asc'
		if (is.null(wt)) { #if wt is null then assume that values in x are the weights
			pos = as.data.frame(which(is.finite(x),arr.ind=TRUE))
			pos$x = getXYcoords(x)$x[pos$row]
			pos$y = getXYcoords(x)$y[pos$col]
			pos$wt = x[cbind(pos$row,pos$col)]
			out = c(COGx=wt.mean(pos$x,pos$wt),COGx.sd=wt.sd(pos$x,pos$wt),COGy=wt.mean(pos$y,pos$wt),COGy.sd=wt.sd(pos$y,pos$wt))
		} else { #if wt is supplied, it must be of the same dim as x and then the values of x are assumed to be your z
			if (!all(dim(x)==dim(wt))) stop('the grids for x & weights must be of the same dimensions')
			pos = as.data.frame(which(is.finite(x),arr.ind=TRUE))
			pos$x = getXYcoords(x)$x[pos$row]
			pos$y = getXYcoords(x)$y[pos$col]
			pos$z = x[cbind(pos$row,pos$col)]
			pos$wt = wt[cbind(pos$row,pos$col)]
			out = c(COGx=wt.mean(pos$x,pos$wt),COGx.sd=wt.sd(pos$x,pos$wt),COGy=wt.mean(pos$y,pos$wt),COGy.sd=wt.sd(pos$y,pos$wt),COGz=wt.mean(pos$z,pos$wt),COGz.sd=wt.sd(pos$z,pos$wt))
		}
	
	}
    # return the output	
	return(out)
}
