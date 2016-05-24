#' @keywords environment
#' @export
#' @title Load a Named Spatial Dataset into the Global Environment
#' @param dataset name of dataset
#' @description Load a named spatial dataset. The package \code{SpatialDataDir} will be
#' searched for the specified dataset.  Set this location with \code{setSpatialDataDir()}.
#' @seealso getSpatialDataDir
#' @seealso setSpatialDataDir
loadSpatialData <- function(dataset=NULL) {
  
  # Use package internal data directory
  dataDir <- getSpatialDataDir()
  
  # Test for and then load the dataset
  filePath <- paste0(dataDir,'/',dataset,'.RData')
  
  if (!file.exists(filePath)) {
    stop(paste0('Spatial data file "',filePath,'" does not exist.',call.=FALSE))
  } else {
    load(filePath,envir=.GlobalEnv)
  }
  
}

