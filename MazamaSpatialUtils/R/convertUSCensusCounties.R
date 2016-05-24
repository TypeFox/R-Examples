#' @keywords datagen
#' @export
#' @title Convert US County Borders Shapefile
#' @param nameOnly logical specifying whether to only return the name without creating the file
#' @description Returns a SpatialPolygonsDataFrame for a US county divisions
#' @details A US county borders shapefile is downloaded and converted to a 
#' SpatialPolygonsDataFrame with additional columns of data. The resulting file will be created
#' in the spatial data directory which is set with \code{setSpatialDataDir()}.
#' @return Name of the dataset being created.
#' @references \url{http://www2.census.gov/geo/tiger/GENZ2013}
#' @seealso setSpatialDataDir
#' @seealso getUSCounty
convertUSCensusCounties <- function(nameOnly=FALSE) {
  
  # Use package internal data directory
  dataDir <- getSpatialDataDir()
  
  # Specify the name of the dataset and file being created
  datasetName <- 'USCensusCounties'
    
  if (nameOnly) return(datasetName)

  # Build appropriate request URL for US County Borders data
  url <- 'http://www2.census.gov/geo/tiger/GENZ2013/cb_2013_us_county_20m.zip'
  
  filePath <- paste(dataDir,basename(url),sep='/')
  utils::download.file(url,filePath)
  # NOTE:  This zip file has no directory so extra subdirectory needs to be created
  utils::unzip(filePath,exdir=paste0(dataDir, '/counties'))
  
  # Convert shapefile into SpatialPolygonsDataFrame
  # NOTE:  The 'counties' directory has been created
  dsnPath <- paste(dataDir,'counties',sep='/')
  shpName <- 'cb_2013_us_county_20m'
  SPDF <- convertLayer(dsn=dsnPath, layerName=shpName, encoding='latin1')
  
  # Rationalize naming:
  # * human readable full nouns with descriptive prefixes
  # * generally lowerCamelCase
  # with internal standards:
  # * countryCode (ISO 3166-1 alpha-2)
  # * stateCode (ISO 3166-2 alpha-2)
  # * longitude (decimal degrees E)
  # * latitude (decimal degrees N)
  # * area (m^2)
  
  # Get STATEFP conversion table from wikipedia. We need this to find state names and codes
  # from STATEFP values.
  # URL of STATEFP conversions
  url <- 'http://en.wikipedia.org/wiki/Federal_Information_Processing_Standard_state_code'

  # Get the raw html from the url
  wikiDoc <- rvest::html(url)
  
  # Get a list of tables in the document
  tables <- rvest::html_nodes(wikiDoc, 'table')
  
  # Assume the relevant list is the first table and parse that into a dataframe
  StateTable <- rvest::html_table(tables[[1]])
  
  # Given a row of country data, use the wikipedia table to find state code and name
  extractState <- function(row, col) {
    state <- StateTable[StateTable['Numeric code']==as.numeric(row['STATEFP']),]
    return(toString(state[col]))
  }
  
  # Standardize naming in the SpatialPolygonsDataFrame
  SPDF$countryCode <- 'US'
  SPDF$countryName <- 'United States'
  SPDF$stateCode <- apply(SPDF@data, 1, extractState, col='Alpha code')
  SPDF$stateName <- apply(SPDF@data, 1, extractState, col='Name')
  SPDF$countyName <- SPDF$NAME
  
  # Subset this dataframe to include only obviously useful columns
  usefulColumns <- c('NAME','ALAND','AWATER','countryCode','countryName','stateCode','stateName',
                     'countyName','COUNTYFP')
  SPDF <- SPDF[,usefulColumns]
  names(SPDF) <- c('name', 'areaLand', 'areaWater', 'countryCode', 'countryName', 
                   'stateCode', 'stateName', 'countyName', 'countyFIPS')
  
  # Group polygons with the same identifier (countyName)
  SPDF <- organizePolygons(SPDF, uniqueID='countyName', sumColumns=c('areaLand','areaWater'))
  
  # Assign a name and save the data
  assign(datasetName,SPDF)
  save(list=c(datasetName),file=paste0(dataDir,'/',datasetName,'.RData'))
  
  # Clean up
  unlink(filePath, force=TRUE)
  unlink(dsnPath, recursive=TRUE, force=TRUE)
  
  return(invisible(datasetName))
}

