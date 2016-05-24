#' Creates *.kml files of a subset of polygons from a polygon shapefile
#' @description Subset polygons acording to the matches between a vector and a specified column from a \link[sp]{SpatialPolygonsDataFrame}.
#' @param shape string with the name of a polygon shapefile or an object of \code{\link{class}} \link[sp]{SpatialPolygonsDataFrame} (see examples).
#' @param psu the values to be matched.
#' @param id column of the *.dbf file with the values to be matched against.
#' @param path string indicating the path to the folder containing the shapfile. If the shapefile is in the working directory or if \code{shape} argument is a shapefile, \code{path} can be ignored.
#' @return *.kml files of the subsetted polygons.
#' @details If there are *.kml files in the working directory, the new created files will overwrite it in case of name matching.
#' 
#' \code{shape} must receive a shapefile with appropriate coordinate reference system, otherwise, \code{MapkmlPSU} report an error.
#' @references \url{http://oswaldosantos.github.io/capm}
#' @seealso \code{\link{readShapeSpatial}}
#' @export
#' @examples
#' # Load data with the polygon identifiers. 
#' data(psu.ssu)
#' 
#' # Take a sample of 10 PSU with probability 
#' # proportional to size with replacement.
#' (selected.psu <- SamplePPS(psu.ssu, 10, write = FALSE))
#' 
#' ## Define shape from shapefile.
#' shp.path <- system.file('extdata', package="capm")
#' # The code above used a shapefile avaliable in the
#' # capm package.
#' # You might want to write a code like:
#' # shp.path  <- 'path_to_the_folder_with_the_shapefile'
#' 
#' # Create *kml files of 10 polygons.
#' # Uncomment the following line to create kml files:
#' # MapkmlPSU('santos', selected.psu[, 1], 1, shp.path)
#' 
#' ## Define the shape argument as an object x of class SpatialPolygonsDataFrame.
#' # MapkmlPSU(x, selected.psu[, 1], 1)
#' 
#' 
MapkmlPSU <- function(shape = NULL, psu = NULL, id = NULL, path = '.') {
  if (class(shape) == 'SpatialPolygonsDataFrame') {
    tmp <- shape
  } else {
    tmp <- readOGR(path, shape)
    #proj4string(tmp) <- CRS('+proj=longlat +ellps=WGS84')
  }
  tmp <- spTransform(tmp, CRS('+proj=longlat +ellps=WGS84'))
  tmp2 = NULL
  for(i in 1:length(psu)) {
    tmp1 <- tmp[which(as.character(tmp@data[ , id]) == psu[i]), ]
    writeOGR(tmp1, dsn = paste(eval(psu[i]), '.kml', sep =''), 
             layer = 'selected_psu', driver = 'KML', overwrite_layer = TRUE)
    tmp2[i] <- which(as.character(tmp@data[ , id]) == psu[i])
  }
  tmp2 <- tmp[tmp2, ]
  if (file.exists('all_psu.kml')) {
    file.remove('all_psu.kml')
  }
  writeOGR(tmp2, dsn = 'all_psu.kml', layer = 'all_selected_psu',
           driver = 'KML', overwrite_layer = TRUE)
  return(cat('\n', 'The maps are in the directory:', '\n\n', getwd()))
}