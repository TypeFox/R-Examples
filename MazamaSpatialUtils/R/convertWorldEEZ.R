#' @keywords datagen
#' @export
#' @title Convert World Exclusive Economic Zones Boundaries Shapefile
#' @param nameOnly logical specifying whether to only return the name without creating the file
#' @description A world EEZ shapefile is downloaded and converted to a 
#' SpatialPolygonsDataFrame with additional columns of data. The resulting file will be created
#' in the spatial data directory which is set with \code{setSpatialDataDir()}.
#' @return Name of the dataset being created.
#' @references \url{http://www.marineregions.org/downloads.php}
#' @seealso setSpatialDataDir
#' @seealso getCountry, getCountryCode
convertWorldEEZ <- function(nameOnly=FALSE) {
  
  # Use package internal data directory
  dataDir <- getSpatialDataDir()
  
  # Specify the name of the file being created
  datasetName <- 'WorldEEZ'
  
  if (nameOnly) return(datasetName)
  
  # Test if the shapefile directory exists.
  if (!file.exists(paste0(dataDir,'/',datasetName))) {
    stop('Shapefile directory does not exists. Please download and convert the shapefile desired.', call.=FALSE)
  } 
  
  # Unzip the downloaded file
  filePath <- paste(dataDir,'World_EEZ_v8_20140228_LR.zip',sep='/')
  utils::unzip(filePath,exdir=paste0(dataDir, '/', datasetName))
  dsnPath <- paste(dataDir, 'WorldEEZ', sep='/')

  # Convert shapefile into SpatialPolygonsDataFrame
  # NOTE:  The 'WorldEEZ' directory has been created
  # NOTE:  Simplify the .shp file using Mapshaper prior to converting layer
  SPDF <- convertLayer(dsn=dsnPath,layerName='World_EEZ_v8_2014')
  
  #   > names(SPDF)
  #   [1] "OBJECTID"    "EEZ"        "Country"    "ID"        "Sovereign"  "Remarks"    "Sov_ID"     "EEZ_ID"     "ISO_3digit"  "MRGID"
  #   [11] "Date_chang" "Area_m2"    "Longitude"  "Latitude"  "MRGID_EEZ"
  
  # Rationalize naming:
  # * human readable full nouns with descriptive prefixes
  # * generally lowerCamelCase
  # with internal standards:
  # * countryCode (ISO 3166-1 alpha-2)
  # * stateCode (ISO 3166-2 alpha-2)
  # * longitude (decimal degrees E)
  # * latitude (decimal degrees N)
  usefulColumns <- c('EEZ', 'Country', 'Sovereign', 'Remarks', 'EEZ_ID', 'ISO_3digit', 'MRGID', 'Area_m2', 'Longitude', 'Latitude')
  SPDF <- SPDF[,usefulColumns]
  names(SPDF) <- c('EEZ', 'country', 'sovereign', 'remarks', 'EEZ_ID', 'ISO3', 'MRGID', 'area', 'longitude', 'latitude')
  
  # Change missing countryCodes to NA
  SPDF$ISO3[SPDF$ISO3 == '-' ] <- NA
  
  # Add more standard columns
  SPDF$countryCode <- codeToCode(SPDF$ISO3)
  SPDF$countryName <- codeToCountry(SPDF$countryCode)
  
  SPDF <- organizePolygons(SPDF, uniqueID='countryCode')
  
  # Assign a name and save the data for World EEZ
  assign(datasetName,SPDF)
  save(list=c(datasetName),file=paste0(dataDir,'/',datasetName,'.RData'))
  
  return(invisible(datasetName))
}

