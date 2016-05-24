#' @title Statistics for edges (lines) based on a defined scale (area).
#' @description Samples rasters for each edge and calculates specified statistics
#'
#' @param x      sp class SpatialLinesDataFrame object 
#' @param r      A rasterLayer, rasterStack or rasterBrick object
#' @param s      Scale (buffer distance) for analysis
#' @param stats  Statistics to calculate. If vectorized, can pass a custom statistic function. 
#' @param echo   Echo progress (TRUE/FALSE)
#'	
#' @return data.frame object with defined statistics, for each raster, at each edge in graph
#'  
#' @examples
#' \dontrun{
#' library(sp)
#' library(spdep)
#' library(raster)
#'   data(rasters)
#'   data(ralu.site)
#'
#' xvars <- stack(rasters[-6])
#' 
#' dist.graph <- knn.graph(ralu.site, row.names = ralu.site@@data[,"SiteName"], 
#'                         max.dist = 500)
#'   str(dist.graph@@data)
#'   
#' skew <- function(x, na.rm = TRUE) {  
#'           if (na.rm) x <- x[!is.na(x)]
#'           sum( (x - mean(x)) ^ 3) / ( length(x) * sd(x) ^ 3 )  
#' 		}
#' 		
#' area.graph.statistics(dist.graph[1:5,], r = xvars, s = 500, 
#'            stats = c("min", "median", "max", "var", "skew")) 
#' }
#' 
#' @export
area.graph.statistics <- function(x, r, s = 500, stats = c("min", "mean", "max"),
                                   echo = FALSE) {
    if(!inherits(x, "SpatialLinesDataFrame")) stop("x is not a SpatialLinesDataFrame object")
	  if (!inherits(r, "RasterLayer") &  
            !inherits(r, "RasterStack") &
		      !inherits(r, "RasterBrick") ) 
	            stop("r is not a raster object")
    bl <- rgeos::gBuffer(x, byid = TRUE, id = rownames(x@data), 
                         width = s, quadsegs = 10)
	results <- as.data.frame(array(0, dim=c(0,length(stats) * dim(r)[3])))
	for( l in 1:length(bl)){
	  if(echo) cat("edge", l, "of", length(bl), "\n")
	  sub.bl <- bl[l,]
	  r.bl <- as.data.frame(raster::extract(r, sub.bl))
	  results.vector <- vector()
	    for(stat in stats) {
          temp.stat <- apply(r.bl, 2, stat)
          names(temp.stat) <- paste(names(r.bl), stat, sep=".")		  
	      results.vector <- append(results.vector, temp.stat) 
	    }
      results <- rbind(results, results.vector)		
	}
	  names(results) <- names(results.vector)
	rownames(results) <- rownames(x@data)
  return(results)	
}
