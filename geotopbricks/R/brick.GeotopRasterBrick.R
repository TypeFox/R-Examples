# TODO: Add comment
# 
# Author: ecor
###############################################################################
NULL

#' 
#' 
#' Added implemetation for 'brick' S4 method 
#' 
#' @param x a 'zoo' object returned by function \code{\link{pointer.to.maps.xyz.time}} or \code{\link{pointer.to.maps.xy.time}} 
#' or a \code{\link{GeotopRasterBrick-class}} object
#' @param layer layer at which raster maps are imported. If is \code{NULL}, maps ara no-zlayer distributed and \code{zoo} must be returend by \code{\link{pointer.to.maps.xy.time}}
#' @param timerange two-elememts vector containing the time range at which geotop maps are imported
#' @param time vector of time instants at which geotop maps are imported
#' @param rows rows of \code{zoo} correspondig to the geotop maps that are imported. By default all rows of \code{zoo} are considered. It is calculated by \code{time} or \code{timerange} if they are not set as \code{NULL}. 
#' @param crs coordinate system see \code{\link{RasterBrick-class}}
#' @param use.read.raster.from.url logical value. Default is \code{TRUE}. If \code{TRUE} the RasterLayer are read with \code{\link{read.raster.from.url}}, istead of \code{\link{raster}} (otherwise). It is recomended in case the files whose paths are contained in \code{x} are remote and are 'http' addresses. In this cases the stand-alone method \code{raster(x)} does not always work and \code{use.read.raster.from.url} is necessary.  
#' @return a  \code{\link{RasterBrick-class}} containing the geopop maps indicated by \code{x}, which is already in a \code{\link{GeotopRasterBrick-class}} object or a 'zoo' object returned by function \code{\link{pointer.to.maps.xyz.time}} or \code{\link{pointer.to.maps.xy.time}}.
#' 
#' @title brick
#' @name brick
#' 
#' @export


#' @import raster
#   @aliases brick
#' @rdname brick-methods
# @keywords methods
# @docType methods
#' @method brick zoo 
#' @aliases brick,zoo-method
#'
#  #' @rdname brick-methods
# @keywords methods
# @docType methods
# @method brick GeotopRasterBrick 
# @aliases brick,GeotopRasterBrick-method
#
#'
#' @examples 
#' # TON TOSS 
#' # See the examples in the functions listed in the 'SeeAlso' section
#' @seealso \code{\link{getvalues.brick.at.depth}},\code{\link{vertical.aggregate.brick.within.depth}}






setMethod('brick', signature(x='zoo'), 
		function(x,layer=1,timerange=NULL,time=NULL,rows=1:nrow(x),crs=NULL,use.read.raster.from.url=TRUE) {
			
			
			if (!is.null(time)) rows <- which(index(x) %in% time)
			if (!is.null(timerange)) rows <- which(index(x)>=timerange[1] & index(x)<=timerange[2])
			
			if ((length(rows)==1) & (length(layer)>1)) {
				 
				 x <- as.matrix(x)
				 x <- x[rows,layer]
			#	 print(rows) getvalues.brick.at.depth raster(x=dtm_map_asc)
			#	 print(layer)
			#	 print(x)
				 list <- as.list(array(NA,length(x)))
				 names(list) <- paste("L",1:length(layer))
				
			} else if (is.null(layer)) {
				
				x <- x[rows]
				list <- as.list(array(NA,length(x)))
				names(list) <- index(x)
				
			} else {
				x <- x[rows,layer]
				list <- as.list(array(NA,length(x)))
				names(list) <- index(x)
			}
			
	
		
			
			if (use.read.raster.from.url) {
				
				for (i in 1:length(x)) {
					
					
					#if (!file.exists(as.character(x[i]))) print(paste("Warning Missing File:",as.character(x[i]),sep=" "))
					# read.raster.from.url	
					list[[i]] <- read.raster.from.url(x=as.character(x[i]))
					
					
				}
				
			} else {
			    for (i in 1:length(x)) {
				
				
				#if (!file.exists(as.character(x[i]))) print(paste("Warning Missing File:",as.character(x[i]),sep=" "))
				# read.raster.from.url	
				list[[i]] <- raster(x=as.character(x[i]))
				
		
			 }
			
			}
			
			
			b <- brick(list)

			if (!is.null(crs)) {
				projection(b) <- crs
			}
			return(b)
		}
)

NULL
#' brick method for GeotopRasterBrick
#' 
#' 
#'  @title brick
#' @name brick
#' 
#' @export
#' @rdname brick-methods
#' @keywords methods
#' @docType methods
#' @method brick GeotopRasterBrick 
#' @aliases brick,GeotopRasterBrick-method

# @export
# @rdname KML-methods
# @keywords methods
# @docType methods
# @method KML GeotopRasterBrick 
# @aliases KML,GeotopRasterBrick-me
setMethod('brick', signature(x='GeotopRasterBrick'), 
		function(x) {
			
			return(x@brick)
		}
)





