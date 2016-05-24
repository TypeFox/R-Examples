#' @title Point sample and statistics for edges (lines)
#' @description Samples rasters for each edge and calculates specified statistics
#'
#' @param x      sp class SpatialLinesDataFrame object 
#' @param r      A rasterLayer, rasterStack or rasterBrick object
#' @param d      Sample distance along edge, for alternate sampline options see sample.line. 
#' @param stats  Statistics to calculate. If vectorized, can pass a custom statistic function. 
#' @param sp     Include an sp class SpatialPointsDataFrame object of the line point samples (FALSE/TRUE) 
#' @param ...    Additional argument passed to sample.line and spsample  
#'	
#' @return data.frame object unless sp=TRUE then, list with data.frame of statistics and SpatialPointsDataFrame contaning point sample data for edges
#'  
#' @examples
#' \dontrun{
#' library(sp)
#' library(spdep)
#' library(raster)
#'   data(rasters)
#'   data(ralu.site)
#'
#' xvars <- stack(rasters)
#' 
#' dist.graph <- knn.graph(ralu.site, row.names = ralu.site@@data[,"SiteName"], 
#'                         max.dist = 1500)
#'   str(dist.graph@data)
#'   
#' skew <- function(x, na.rm = TRUE) {  
#'           if (na.rm) x <- x[!is.na(x)]
#'           sum( (x - mean(x)) ^ 3) / ( length(x) * sd(x) ^ 3 )  
#' 		}
#' 		
#' stats <- graph.statistics(dist.graph, r = xvars, d=30, 
#'             stats = c("min", "median", "max", "var", "skew"),
#'             sp = FALSE) 
#' 			
#' dist.graph@@data <- data.frame(dist.graph@@data, stats)
#'   str(dist.graph@@data)
#' }
#' 
#' @export
graph.statistics <- function(x, r, d = 30, stats = c("min", "mean", "max"), 
                             sp = FALSE, ...) {
    if(!inherits(x, "SpatialLinesDataFrame")) stop("x is not a SpatialLinesDataFrame object")
	  if (!inherits(r, "RasterLayer") &  
            !inherits(r, "RasterStack") &
		      !inherits(r, "RasterBrick") ) 
	            stop("r is not a raster object")
    dots <- as.list(match.call(expand.dots = TRUE)[-1])
    dots[["x"]] <- x
  if (is.null(dots[["min.samp"]]) & "min.samp" %in% names(dots) == FALSE) dots[["min.samp"]] <-  2
  if (is.null(dots[["type"]]) & "type" %in% names(dots) == FALSE) dots[["type"]] <-  "regular"
  samp <- do.call(spatialEco::sample.line, dots)
  samp@data <- data.frame(samp@data, raster::extract(r, samp))
  results <- NULL
    if(dim(r)[3] > 1) {   
      for( p in names(r) ) {
        for(s in 1:length(stats)) {
          results <- cbind(results, as.vector(tapply(samp@data[,p], samp@data[,"LID"], stats[s])))
        }
	  }
	} else {
	  names(samp@data)[2] <- names(r)
      for(s in 1:length(stats)) {
          results <- cbind(results, as.vector(tapply(samp@data[,2], samp@data[,"LID"], stats[s])))
      } 
    }	  
    results <- as.data.frame(results)
	  rn <- vector()
	    for(n in names(r)) { rn <- append(rn, paste(stats, n, sep="."))}
	      names(results) <- rn
  if( sp == TRUE) {
    return(list(statistics = results, sample = samp) )
  } else {
    return( results )
  }  
}
