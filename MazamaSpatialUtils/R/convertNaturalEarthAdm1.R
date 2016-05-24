#' @keywords datagen
#' @export
#' @title Convert Level 1 (State) Borders Shapefile
#' @param nameOnly logical specifying whether to only return the name without creating the file
#' @description Returns a SpatialPolygonsDataFrame for a 1st level administrative divisions
#' @details A state border shapefile is downloaded and converted to a 
#' SpatialPolygonsDataFrame with additional columns of data. The resulting file will be created
#' in the spatial data directory which is set with \code{setSpatialDataDir()}.
#' 
#' Within the \pkg{MazamaSpatialUtils} package the phrase 'state' refers to administrative
#' divisions beneath the level of the country or nation. This makes sense in the United 'States'. In
#' other countries this level is known as 'province', 'territory' or some other term.
#' @return Name of the dataset being created.
#' @references \url{http//www.naturalearthdata.com/download}
#' @references \url{http://www.statoids.com/ihasc.html}
#' @seealso setSpatialDataDir
#' @seealso getState, getStateCode
convertNaturalEarthAdm1 <- function(nameOnly=FALSE) {
  
  # Use package internal data directory
  dataDir <- getSpatialDataDir()
  
  # Specify the name of the dataset and file being created
  datasetName <- 'NaturalEarthAdm1'
    
  if (nameOnly) return(datasetName)
  
  # Specify administrative levels for URL
  adm <- 1
  level <- 'states_provinces'

  # Build appropriate request URL for Natural Earth data
  url <- paste0('http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_',
                adm, '_',
                level, '.zip')

  filePath <- paste(dataDir,basename(url),sep='/')
  utils::download.file(url,filePath)
  # NOTE:  This zip file has no directory so extra subdirectory needs to be created
  utils::unzip(filePath,exdir=paste0(dataDir, '/adm'))
    
  # Convert shapefile into SpatialPolygonsDataFrame
  # NOTE:  The 'adm' directory has been created
  dsnPath <- paste(dataDir,'adm',sep='/')
  shpName <- paste('ne', '10m_admin', adm, level, sep='_')
  SPDF <- convertLayer(dsn=dsnPath, layerName=shpName)

  # Rationalize naming:
  # * human readable full nouns with descriptive prefixes
  # * generally lowerCamelCase
  # with internal standards:
  # * countryCode (ISO 3166-1 alpha-2)
  # * stateCode (ISO 3166-2 alpha-2)
  # * longitude (decimal degrees E)
  # * latitude (decimal degrees N)
  # * area (m^2)
  
  # Add the core identifiers to the SpatialPolygonsDataFrame
  SPDF$countryCode <- SPDF$iso_a2
  SPDF$stateCode <- stringr::str_split_fixed(SPDF$code_hasc,'\\.',5)[,2]
  SPDF$countryName <- MazamaSpatialUtils::codeToCountry(SPDF$countryCode)
  SPDF$stateName <- SPDF$name
  
  # Subset this dataframe to include only obviously useful columns
  usefulColumns <- c('countryCode','countryName','stateCode','stateName','latitude','longitude','area_sqkm',
                     'postal','adm1_code','code_hasc','fips','gns_lang','gns_adm1')
  # TODO:  Check that usefulColumns are actually found in the dataframe
  SPDF <- SPDF[,usefulColumns]
  
  # Rationalize units:
  # * SI  
  # NOTE:  Area seems to be in units of km^2. Convert these to m^2
  SPDF$area <- SPDF$area_sqkm * 1e6  

  # Use NA instead of -99 and -90 for missing values
  SPDF@data[SPDF@data == -99] <- NA
  SPDF@data[SPDF@data == -90] <- NA
  
  # Group polygons with the same identifier (gns_adm1)
  SPDF <- organizePolygons(SPDF, uniqueID='adm1_code', sumColumns=c('area','area_sqkm'))  
  
  # Assign a name and save the data
  assign(datasetName,SPDF)
  save(list=c(datasetName),file=paste0(dataDir,"/",datasetName,'.RData'))
  
  # Clean up
  unlink(filePath, force=TRUE)
  unlink(dsnPath, recursive=TRUE, force=TRUE)
  
  return(invisible(datasetName))
}

