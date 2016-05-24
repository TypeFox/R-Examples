#' Get a catchment object from selected shape file.
#' @param filePath A string representing the path of the shape file.
#' @return A catchment object can be used in \code{getSpatialMap()}.
#' @export
#' @details This function is based on the package \code{rgdal} and \code{sp}, and the output comes from the package 
#' \code{sp}
#' @examples
#' #open internal file
#' file <- system.file("extdata", "testCat.shp", package = "hyfo")
#' catchment <- shp2cat(file)
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @import rgdal
#' @importFrom utils tail
#' @references 
#' 
#' \itemize{
#' \item Roger Bivand, Tim Keitt and Barry Rowlingson (2015). rgdal: Bindings for the Geospatial Data
#' Abstraction Library. R package version 1.0-4. http://CRAN.R-project.org/package=rgdal
#' 
#' \item R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' }
#' 
#' 
shp2cat <- function(filePath) {
  #if the path <- file.choose(), the seperator is '\\'
  if (grepl('\\\\', filePath)) {
    catName <- tail(strsplit(filePath,'\\\\')[[1]], 1)#needs to be four \, caused by some window system problem
    catName1 <- strsplit(catName, '\\.')[[1]][1]
    catName2 <- paste('\\\\', catName, sep = '')
    folderName <- strsplit(filePath, catName2)[[1]]
    n <- list.files(folderName, pattern = catName1)
    if (length(n) == 1) stop('Please place the shp file in the folder containing 
                             full related files, not only the shape file')
  #the other seperator is '/'  
  } else if (grepl('/', filePath)) {
    catName <- tail(strsplit(filePath,'/')[[1]], 1)#needs to be four \, caused by some window system problem
    catName1 <- strsplit(catName, '\\.')[[1]][1]
    catName2 <- paste('/', catName, sep = '')
    folderName <- strsplit(filePath, catName2)[[1]]
    n <- list.files(folderName, pattern = catName1)
    if (length(n) == 1) stop('Please place the shp file in the folder containing 
                             full related files, not only the shape file')
  }
  
  if (length(folderName) == 0) stop('No shape file found, make sure the shp file is selected.')
  catchment <- readOGR(folderName, catName1)
  return(catchment)
}
