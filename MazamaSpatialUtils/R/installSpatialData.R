#' @keywords environment
#' @export
#' @title Install a Named Spatial Dataset
#' @param dataset name of dataset
#' @param verbose logical flag controlling detailed progress statements
#' @param ... additional arguments needed by some \code{convert~} functions
#' @description Install a named spatial dataset.
#' @return Character name of the installed dataset.
installSpatialData <- function(dataset=NULL, verbose=FALSE, ...) {
  
  # Use package internal data directory
  dataDir <- getSpatialDataDir()
  
  # Get the appropriate 'convert' function
  FUN <- get(paste0('convert',dataset))

  # Determine the file name and absolute path
  datasetName <- FUN(..., nameOnly=TRUE)
  filePath <- paste0(dataDir,'/',datasetName,'.RData')
  
  if (file.exists(filePath)) {
    if (verbose) message(paste0(filePath,' already exists.'))
  } else {
    # Download/Convert/Install the dataset
    FUN(..., nameOnly=FALSE)    
  }
  
}

