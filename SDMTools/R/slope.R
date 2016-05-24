#' Slope and aspect calculations
#' 
#' \code{slope} and \code{aspect} calculates the slope and aspect of raster
#' surfaces of class 'asc' (SDMTools & adehabitat packages), 'RasterLayer'
#' (raster package) or 'SpatialGridDataFrame' (sp package).\cr\cr Methods are
#' based on Burrough and McDonell (1998).
#' 
#' Slope returns values representing the 'rise over run' with "run" units
#' representing cellsize if \code{latlon}=FALSE or km if \code{latlon}=TRUE.
#' This can be changed to percentage (multiply by 100) or to degrees by ATAN (
#' \code{output} ) * 57.29578.\cr\cr Aspect returns the direction (0 to 360)
#' with North being 0. Values of -1 are flat areas with no slope or
#' aspect.\cr\cr As this method requires information from the surrounding
#' cells, missing data (NAs or edges) are populated with the value from the
#' 'cell-of-interest').
#' 
#' @param mat a matrix of data representing z heights. Matrix can be a raster
#' of class 'asc' (adehabitat package), 'RasterLayer' (raster package) or
#' 'SpatialGridDataFrame' (sp package)
#' @param latlon boolean value representing if the data is geographic.
#' @return an object of the same class as \code{mat}.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @references Burrough, P. A. and McDonell, R.A., 1998. Principles of
#' Geographical Information Systems (Oxford University Press, New York), p.
#' 190.
#' @examples
#' #define a simple asc with some slope and direction
#' tasc = as.asc(matrix(1:50,nr=10,nc=5),yll=75); tasc[,]
#' slope(tasc)[,] #show the output of slope
#' aspect(tasc)[,] #show the output of the aspect
#' 
#' #define a FLAT simple asc 
#' tasc = as.asc(matrix(10,nr=10,nc=5),yll=75); tasc[,]
#' slope(tasc)[,] #show the output of slope
#' aspect(tasc)[,] #show the output of the aspect
#' 
#' 
#' @export 
#' @useDynLib SDMTools Slope Aspect
slope <- 
function(mat,latlon=FALSE) {
	#check input for class for returning info
	if (any(class(mat) == 'asc')) { attrib = attributes(mat)
	} else if (any(class(mat) %in% 'RasterLayer')) { attrib = mat; mat = asc.from.raster(mat)
	} else if (any(class(mat) == 'SpatialGridDataFrame')) { attrib = mat; mat = asc.from.sp(mat)
	} else { attrib = attributes(mat) }
	if (!any(class(mat) == 'asc')) { stop('objects must be of class "asc"') } #check to ensure asc 
	
	# get the cell size information
	if (latlon) {
		tt = grid.info(getXYcoords(mat)$y,attr(mat,'cellsize')) #if latlon = true get the length & width of cells
		width = rowMeans(cbind(tt$top,tt$bottom))/1000; height = tt$side/1000 #get the width and height of the cells in km (NOT m)
	} else { width = height = rep(attr(mat,'cellsize'),length(getXYcoords(mat)$y)) } #get the cell width & height 
	
	slop = t(mat[,dim(mat)[2]:1])
	slop = .Call('Slope',slop,width,height,PACKAGE='SDMTools') #get the slope information
	mat[,] = t(slop[dim(slop)[1]:1,]) #move all slope info to mat
	
	#reset the attributes of the input
	if (any(class(attrib) %in% 'RasterLayer')) { attrib = setValues(attrib, as.vector(t(t(unclass(mat))[dim(mat)[2]:1,]))); return(attrib)
	} else if (any(class(attrib) == 'SpatialGridDataFrame')) { attrib@data[1] = as.vector(unclass(mat)[,dim(mat)[2]:1]); return(attrib)
	} else { attributes(mat) = attrib; return(mat) }
}

#' @rdname slope
#' @export
aspect <- 
function(mat,latlon=FALSE) {
	#check input for class for returning info
	if (any(class(mat) == 'asc')) { attrib = attributes(mat)
	} else if (any(class(mat) %in% 'RasterLayer')) { attrib = mat; mat = asc.from.raster(mat)
	} else if (any(class(mat) == 'SpatialGridDataFrame')) { attrib = mat; mat = asc.from.sp(mat)
	} else { attrib = attributes(mat) }
	if (!any(class(mat) == 'asc')) { stop('objects must be of class "asc"') } #check to ensure asc 
	
	# get the cell size information
	if (latlon) {
		tt = grid.info(getXYcoords(mat)$y,attr(mat,'cellsize')) #if latlon = true get the length & width of cells
		width = rowMeans(cbind(tt$top,tt$bottom))/1000; height = tt$side/1000 #get the width and height of the cells in km (NOT m)
	} else { width = height = rep(attr(mat,'cellsize'),length(getXYcoords(mat)$y)) } #get the cell width & height 
	
	asp = t(mat[,dim(mat)[2]:1]) #reset the grid so that [1,1] is the North West corner (not the default of lower-left with reversed lat & lon)
	asp = .Call('Aspect',asp,width,height,PACKAGE='SDMTools') #get the aspect information
	mat[,] = t(asp[dim(asp)[1]:1,]) #move all aspect info to mat

	#reset the attributes of the input
	if (any(class(attrib) %in% 'RasterLayer')) { attrib = setValues(attrib, as.vector(t(t(unclass(mat))[dim(mat)[2]:1,]))); return(attrib)
	} else if (any(class(attrib) == 'SpatialGridDataFrame')) { attrib@data[1] = as.vector(unclass(mat)[,dim(mat)[2]:1]); return(attrib)
	} else { attributes(mat) = attrib; return(mat) }
}
