#' Landscape Zonal Statistics
#' 
#' \code{ZonalStat} calculates the statistics of data for specified zones of
#' two matrices of data. The matrix can be a raster of class 'asc' (adehabitat
#' package), 'RasterLayer' (raster package) or 'SpatialGridDataFrame' (sp
#' package).
#' 
#' The code summarizes the data for defined zones. Nearly any function can be
#' used for summarizing the data. \cr \cr The FUN defined with 'all' as one of
#' or the only function will append the functions of min, 1st quarter, median,
#' 3rd quarter, max, mean, standard deviation and n to what is being
#' calculated.
#' 
#' @param mat a matrix of data to be summarized; The matrix can be a raster of
#' class 'asc' (adehabitat package), 'RasterLayer' (raster package) or
#' 'SpatialGridDataFrame' (sp package)
#' @param zones a matrix of data with individual patches identified as with
#' \code{ConnCompLabel}; The matrix must be of the same size & extent as
#' \code{mat}
#' @param FUN a single or vector of functions to be applied to each 'zone'; the
#' default of 'all' will calculate min, 1st quarter, median, 3rd quarter, max,
#' mean, standard deviation and n
#' @return a data.frame listing \item{zone}{the unique ID for each zone.}
#' \item{functions...}{a column for each of the functions identified}
#' 
#' The data.frame will have an atribute defining the number of NA values that
#' were excluded from the analysis.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @examples
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
#' ccl.mat #this is the zone matrix to be used
#' 
#' #create a random data matrix
#' data.mat = matrix(runif(100),nr=10,nc=10)
#' data.mat
#' 
#' #calculate the zonal statistics
#' zs.data = ZonalStat(data.mat,ccl.mat,FUN='all')
#' zs.data
#' 
#' #just calculate the sum
#' zs.data = ZonalStat(data.mat,ccl.mat,FUN='sum')
#' zs.data
#' 
#' #calculate sum & n & 'all' and show when a function is not defined
#' zs.data = ZonalStat(data.mat,ccl.mat,
#'     FUN=c('sum','length','not.a.function','all'))
#' zs.data
#' attr(zs.data,'excluded NAs') #show how many NAs were omitted from analysis
#' 
#' @export
ZonalStat <- 
function(mat,zones,FUN='all')	{
	#check if functions are defined
	for (fun in FUN) { 
		if (!is.function(try(match.fun(fun),silent=TRUE))) { 
			FUN = FUN[-which(FUN==fun)]; warning(paste(fun,' is not a defined function!'))
		}
	}
	if (length(FUN)<1) stop('no known functions defined in request')

	#check if raster from sp or raster package and convert if necessary
	if (any(class(mat) %in% 'RasterLayer')) mat = asc.from.raster(mat)
	if (any(class(mat) == 'SpatialGridDataFrame')) mat = asc.from.sp(mat)
	if (any(class(zones) %in% 'RasterLayer')) zones = asc.from.raster(zones)
	if (any(class(zones) == 'SpatialGridDataFrame')) zones = asc.from.sp(zones)
	#check to ensure matrix
	mat = try(as.matrix(mat)); if (!is.matrix(mat)) stop('objects must be a matrix')
	zones = try(as.matrix(zones)); if (!is.matrix(zones)) stop('objects must be a matrix')
	if (!(dim(zones)[1]==dim(mat)[1] & dim(zones)[2]==dim(mat)[2])) stop('objects must be of the same extent / cellsize')
	#setup a matrix to summarize
	tt = data.frame(zones=as.vector(zones),data=as.vector(mat)); na.count = nrow(tt)
	tt = na.omit(tt); na.count = na.count-nrow(tt) #omit and count NA values
	tt = split(tt$data,tt$zones)
	
	#cycle through each of the functions and extract the data
	if (any(FUN=='all')) { FUN = c(FUN,'quantile','mean','sd','length'); FUN=FUN[-which(FUN=='all')]; FUN=unique(FUN) }
	
	#define the output
	out = data.frame(zones=names(tt))
	for (fun in FUN) { 
		if (fun == 'quantile') {
			out$min = out$qtr.25 = out$median = out$qtr.75 = out$max = NA
		} else { out[[fun]] = NA }		
	}
	
	#cycle through each of the groups in tt
	for (ii in 1:length(tt)) {
		#cycle through each of the functions
		for (fun in FUN) {
			if (fun=='quantile') {
				qq = quantile(tt[[ii]])
				out$min[ii] = qq[1]; out$qtr.25[ii] = qq[2]; out$median[ii] = qq[3]; out$qtr.75[ii] = qq[4]; out$max[ii] = qq[5]
			} else { out[ii,fun] = match.fun(fun)(tt[[ii]]) }
		}
	}
	#add an attribute defining how many NAs were removed
	attr(out,'excluded NAs') = na.count
	
	#return the data
	return(out)
}
