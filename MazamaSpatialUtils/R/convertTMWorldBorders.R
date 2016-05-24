#' @keywords datagen
#' @export
#' @title Convert World Borders Shapefile
#' @param nameOnly logical specifying whether to only return the name without creating the file
#' @description Returns a SpatialPolygonsDataFrame for world divisions
#' @details A world borders shapefile is downloaded and converted to a 
#' SpatialPolygonsDataFrame with additional columns of data. The resulting file will be created
#' in the spatial data directory which is set with \code{setSpatialDataDir()}.
#' @return Name of the dataset being created.
#' @references \url{http://thematicmapping.org/downloads/}
#' @seealso setSpatialDataDir
#' @seealso getCountry, getCountryCode
convertTMWorldBorders <- function(nameOnly=FALSE) {

  # Use package internal data directory
  dataDir <- getSpatialDataDir()
    
  # Specify the name of the dataset and file being created
  datasetName <- 'TMWorldBorders'
  
  if (nameOnly) return(datasetName)
  
  # Build appropriate request URL for TM World Borders data
  url <- 'http://thematicmapping.org/downloads/TM_WORLD_BORDERS-0.3.zip'
  
  filePath <- paste(dataDir,basename(url),sep='/')
  utils::download.file(url,filePath)
  # NOTE:  This zip file has no directory so extra subdirectory needs to be created
  utils::unzip(filePath,exdir=paste0(dataDir,'/world'))
  
  # Convert shapefile into SpatialPolygonsDataFrame
  # NOTE:  The 'world' directory has been created
  dsnPath <- paste(dataDir,'world',sep='/')
  shpName <- 'TM_WORLD_BORDERS-0.3'
  SPDF <- convertLayer(dsn=dsnPath,layerName=shpName)

  # Rationalize naming:
  # * human readable full nouns with descriptive prefixes
  # * generally lowerCamelCase
  # with internal standards:
  # * countryCode (ISO 3166-1 alpha-2)
  # * stateCode (ISO 3166-2 alpha-2)
  # * longitude (decimal degrees E)
  # * latitude (decimal degrees N)
  # * area (m^2)
  
  # Relabel and standardize the naming in the SpatialPolygonsDataFrame
  names(SPDF) <- c('FIPS','countryCode','ISO3','UN_country','countryName','area','population2005','UN_region','UN_subregion','longitude','latitude')
  
  # Rationalize units:
  # * SI  
  # NOTE:  Area seems to be in units of (10 km^2). Convert these to m^2
  SPDF$area <- SPDF$area * 1e7
  
  # Group polygons with the same identifier (countryCode)
  SPDF <- organizePolygons(SPDF, uniqueID='countryCode', sumColumns=c('area','population2005'))
  # NOTE:  This dataset already has grouped polygons
  
  # Assign a name and save the data
  assign(datasetName,SPDF)
  save(list=c(datasetName),file=paste0(dataDir,'/',datasetName,'.RData'))
  
  # Clean up
  unlink(filePath, force=TRUE)
  unlink(dsnPath, recursive=TRUE, force=TRUE)
  
  return(invisible(datasetName))
}

