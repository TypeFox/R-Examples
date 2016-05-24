# # TODO: Add comment
# # 
# # Author: johnsond@afsc.noaa.gov
# ###############################################################################
# 
# 
# 
# #' Compute a spatial use grid from a crawl prediction
# #' 
# #' This function take a SpatialPoints object and a spatial GridTopology
# #' object from the 'sp' package and outputs the number of point locations
# #' in each grid cell
# #' 
# #' 
# #' @param object A 'SpatialPoints' object created with the \code{sp} package.  T
# #' @param grid A \code{GridTopology} object from the 'sp' package
# #' @param subset An indicator of which times should be used for calculation of
# #' the use grid. Can be a logical vector or a vector of integers, such as from
# #' a call to \code{which}
# #' @return A \code{SpatialGridDataFrame} with data column 'use' which gives the
# #' count of locations within the grid cell
# #' @author Devin S. Johnson <devin.johnson@@noaa.gov>
# #' @export
# #' @import sp 
# #' @import raster 
# crwUseGrid <- function(object, grid, subset=TRUE){
#   object <- object[subset,]
#   useTemplate <- SpatialGrid(grid=grid, proj4string=CRS(proj4string(object)))
#   out <- as(rasterize(object, raster(useTemplate), fun=sum), "SpatialGridDataFrame")
#   names(out@data) <- "use"
#   return(out)
# }  
#   