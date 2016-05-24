#' Create a Grid of Cell Areas or Perimeters
#' 
#' Creates a grid of cell areas or perimeters for spatial grids in geographic
#' (lat-lon) projections.
#' 
#' 
#' @param mat a matrix representing a raster of class 'asc' (this & adehabitat
#' package), 'RasterLayer' (raster package) or 'SpatialGridDataFrame' (sp
#' package)
#' @return \item{grid.area}{Returns an ascii grid file which contains the
#' values of the area in each cell.} \item{grid.perimter}{Returns an ascii grid
#' file which contains the values of the perimeter in each cell. }
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com} & Lorena Falconi
#' \email{lorefalconi@@gmail.com}
#' @examples
#' 
#' #Create an ascii file
#' y=seq(10,50,0.5)
#' x=seq(140,180,0.5)
#' cellsize=0.5
#' data1=sample(160,140)
#' out1.asc=as.asc(matrix(data1,nc=y, nr=x), xll=min(x), yll=min(y), cellsize=cellsize)
#' 
#' grid.area(out1.asc)[,]
#' 
#' grid.perimeter(out1.asc)[,]
#' 
#' @export 
grid.area <-
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
	#check to ensure asc 
	if (!any(class(mat) == 'asc')) { stop('objects must be of class "asc"') }
	#apply the gridinfo 
	area = grid.info(getXYcoords(mat)$y,attr(mat,'cellsize'))$area
	mat[is.finite(mat)] = 1; for (ii in 1:length(area)) mat[,ii] = mat[,ii] * area[ii]
	#reset the attributes of the input
	if (any(class(attrib) %in% 'RasterLayer')) {
		attrib = setValues(attrib, as.vector(t(t(unclass(mat))[dim(mat)[2]:1,]))); return(attrib)
	} else if (any(class(attrib) == 'SpatialGridDataFrame')) {
		attrib@data[1] = as.vector(unclass(mat)[,dim(mat)[2]:1]); return(attrib)
	} else {
		attributes(mat) = attrib; return(mat)
	}
}

#' @rdname grid.area
#' @export
grid.perimeter <-
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
	#check to ensure asc 
	if (!any(class(mat) == 'asc')) { stop('objects must be of class "asc"') }
	#apply the gridinfo 
	perim = grid.info(getXYcoords(mat)$y,attr(mat,'cellsize'))
	perim = perim$top+perim$bottom+2*perim$side
	mat[is.finite(mat)] = 1; for (ii in 1:length(perim)) mat[,ii] = mat[,ii] * perim[ii]
	#reset the attributes of the input
	if (any(class(attrib) %in% 'RasterLayer')) {
		attrib = setValues(attrib, as.vector(t(t(unclass(mat))[dim(mat)[2]:1,]))); return(attrib)
	} else if (any(class(attrib) == 'SpatialGridDataFrame')) {
		attrib@data[1] = as.vector(unclass(mat)[,dim(mat)[2]:1]); return(attrib)
	} else {
		attributes(mat) = attrib; return(mat)
	}
}
