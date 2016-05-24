#' Least Cost Moving Windows Calculation
#' 
#' This is a moving window that for each cell returns the minimum 'cost' based
#' on surrounding data cells and some dispersal distance cost.
#' 
#' This method moves over the matrix of values, summing the moving window cost
#' \code{mw} and the matrix \code{mat}, returning the minimum cost value. This
#' was created to estimate the least cost path through time for all cells in a
#' matrix (see example).
#' 
#' @param mat a matrix of values that can be based on a raster dataset. Lower
#' values should represent lower cost. The matrix can be a raster of class
#' 'asc' (adehabitat package), 'RasterLayer' (raster package) or
#' 'SpatialGridDataFrame' (sp package)
#' @param mw a distance-cost matrix to be applied to each cell of 'mat'. This
#' matrix can be dispersal costs. Lower values should represent lower cost.
#' @param mnc an integer value representing the radius for 'mw' in number of
#' cells.
#' @return A matrix of values of the same dimensions and class as input
#' \code{mat}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @examples
#' 
#' 
#' #create a simple object of class 'asc'
#' tasc = as.asc(matrix(1:100,nr=10,nc=10)); print(tasc)
#' 
#' #show the input matrix
#' print(tasc[1:10,1:10])
#' 
#' #vary the moving windows
#' 
#' ###no cost window of 2 cell radius
#' tcost = matrix(0,nr=5,nc=5); print(tcost)
#' out = lcmw(tasc, tcost, 2); print(out[1:10,1:10])
#' 
#' ###no cost with a circular radius of 2
#' tcost = matrix(NA,nr=5,nc=5)
#' #populate the distances
#' for (y in 1:5){
#'     for (x in 1:5){
#'         tcost[y,x] = sqrt((3-y)^2 + (3-x)^2)
#'     }
#' }
#' 
#' #remove distance values > max.num.cells
#' tcost[which(tcost>2)]=NA
#' 
#' #no cost matrix
#' tcost1 = tcost; tcost1[is.finite(tcost1)]=1; print(tcost1)
#' out = lcmw(tasc, tcost1, 2); print(out[1:10,1:10])
#' 
#' #linear cost
#' tcost = tcost/2; print(tcost)
#' out = lcmw(tasc, tcost, 2); print(out[1:10,1:10])
#' 
#' 
#' @export 
#' @useDynLib SDMTools getmin movewindow
lcmw <-
function(mat,mw,mnc) {
	#check input for class for returning info
	if (class(mat) == 'asc') { 
		attrib = attributes(mat)
	} else if (any(class(mat) %in% 'RasterLayer')) {
		attrib = mat; mat = asc.from.raster(mat)
	} else if (any(class(mat) == 'SpatialGridDataFrame')) {
		attrib = mat; mat = asc.from.sp(mat)
	} else {
		attrib = attributes(mat)
	}
	#buffer edges by full number of distance cells
	#define the shifts in mat to account for a moving window...
	vals = expand.grid(Y=-mnc:mnc,X=-mnc:mnc) #define all shifts
	vals$cost = mw[(mnc+1)+cbind(vals$Y,vals$X)];vals=na.omit(vals) #extract the cost of associated with the move
	nrow.vals = nrow(vals)
	#cycle through and get the output
	if (nrow.vals <5000) {
		return(.Call("movewindow",mat,as.integer(vals$X),as.integer(vals$Y),as.numeric(vals$cost),PACKAGE='SDMTools'))
	} else {
		num.subsets = nrow.vals%/%2000
		#run the first set of 2000
		tmin = 1; tmax = 2000
		#print a status
		cat('0%...')
		#create the first part of the moving window
		out = .Call("movewindow",mat,as.integer(vals$X[tmin:tmax]),as.integer(vals$Y[tmin:tmax]),as.numeric(vals$cost[tmin:tmax]),PACKAGE='SDMTools')
		#cycle through the remaining data
		for (i in 1:num.subsets){
			if (i<num.subsets){
				tmin = i*2000+1; tmax = (i+1)*2000
			} else {
				tmin = i*2000+1; tmax = nrow.vals
			}
			cat(round(tmin/nrow.vals*100,1),'%...',sep='')
			out2 = .Call("movewindow",mat,as.integer(vals$X[tmin:tmax]),as.integer(vals$Y[tmin:tmax]),as.numeric(vals$cost[tmin:tmax]),PACKAGE='SDMTools')
			out = .Call("getmin",out,out2,PACKAGE='SDMTools')
			if (dim(out)[1] != dim(mat)[1] | dim(out)[2] != dim(mat)[2]) print('error in dimensions...check output')
		}
		cat('done\n')
	}
	#reset the attributes of the input
	if (any(class(attrib) %in% 'RasterLayer')) {
		attrib = setValues(attrib, as.vector(t(t(unclass(out))[dim(out)[2]:1,]))); return(attrib)
	} else if (any(class(attrib) == 'SpatialGridDataFrame')) {
		attrib@data[1] = as.vector(unclass(out)[,dim(out)[2]:1]); return(attrib)
	} else {
		attributes(out) = attrib; return(out)
	}
}

