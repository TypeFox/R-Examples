#' @keywords datagen
#' @export
#' @title Convert USGS Hydrologic Unit Shapefiles
#' @param dsnPath directory where the WBD HUC datasets are found
#' @param level character or integer which must be 2, 4, 6, 8, 10, 12 or 14
#' @param extension character extsion associated with mapshaper simplified files
#' @param nameOnly logical specifying whether to only return the name without creating the file
#' @description Previously downloaded shapefiles from the USGS 
#' \href{http://nhd.usgs.gov/wbd.html}{Watershed Boundary Dataset} are converted to a 
#' SpatialPolygonsDataFrame with additional columns of data. The resulting file will be
#' created in the spatial data directory which is set with \code{setSpatialDataDir()}.
#' @details The full WBD dataset can be downloaded from the USGS with the 
#' following command:
#' \preformatted{
#' curl ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/WBD/Shape/WBD_National.zip -O
#' }
#' 
#' Typically, the raw data will be simplified using the command line version of
#' \href{https://github.com/mbloch/mapshaper}{mapshpaper}. (Installation instructions are
#' found at this URL.)
#' 
#' With mapshaper, you can reduce the number of vertices in the polygons, greatlyl improving
#' the efficiency of spatial searches. Experimentation at the 
#' \href{http://www.mapshaper.org}{mapshaper website} show that a reduction to 1-2%
#' of the original shapefile size still retains the recognizable shape of polygons, removing
#' only the higher order "crenellations" in the polygons.
#' 
#' An example use of mapshaper would be:
#' \preformatted{
#' mapshaper WBDHU2.shp --simplify 1% --o WBDHU2_01.shp
#' }
#' 
#' A full suite of \code{.shp, .shx, .dbf, .prj} files will be created for the new name \code{WBDHU2_02}.
#' 
#' @return Name of the dataset being created.
#' @references \url{http://nhd.usgs.gov/wbd.html}
#' @seealso setSpatialDataDir

# TODO:  Convert missing state codes to state codes from allStateCode, with a note explaining 
# TODO:  how and why. Figure out why it is printing all those numbers when it runs and change.  

convertWBDHUC <- function(dsnPath=NULL, level=8, extension="", nameOnly=FALSE) {
  
  # Sanity check dsnPath
  if ( is.null(dsnPath) ) stop(paste0('Argument dsnPath must be specified.'))
  if ( !file.exists(dsnPath) ) stop(paste0('dsnPath="',dsnPath,'" not found.'))
  
  # 'level' should be a character string
  level <- as.character(level)
  
  # Use package internal data directory
  dataDir <- getSpatialDataDir()
  
  # Specify the name of the dataset and file being created
  datasetName <- paste0('WBDHU', level) 
  
  if (nameOnly) return(datasetName)

  # Convert shapefile into SpatialPolygonsDataFrame
  layerName <- paste0('WBDHU', level, extension)
  SPDF <- convertLayer(dsn=dsnPath, layerName=layerName)

  # Rationalize naming:
  # * human readable full nouns with descriptive prefixes
  # * generally lowerCamelCase
  # with internal standards:
  # * countryCode (ISO 3166-1 alpha-2)
  # * stateCode (ISO 3166-2 alpha-2)
  # * longitude (decimal degrees E)
  # * latitude (decimal degrees N)
  
  # Subset this dataframe to include only obviously useful columns
  
  # NOTE:  Comments are relevant to the WBD as downlaoded on 2015-12-04
  
  if ( level == '2' ) {
    # 22 features
    #
    # [1] "SHAPE_AREA" "SOURCEFEAT" "AREASQKM"   "METASOURCE" "SOURCEORIG" "LOADDATE"   "TNMID"     
    # [8] "AREAACRES"  "GNIS_ID"    "NAME"       "SOURCEDATA" "STATES"     "HUC2"       "SHAPE_LENG"
    #[15] "GEODB_OID"  "OBJECTID"  
    usefulColumns <- c('AREASQKM', 'HUC2', 'NAME', 'STATES')
  } else if ( level == '4' ) {
    # 223 features
    #
    # [1] "AREAACRES"  "AREASQKM"   "GNIS_ID"    "HUC4"       "LOADDATE"   "METASOURCE" "NAME"      
    # [8] "SHAPE_AREA" "SHAPE_LENG" "SOURCEDATA" "SOURCEFEAT" "SOURCEORIG" "STATES"     "TNMID"     
    usefulColumns <- c('AREASQKM', 'HUC4', 'NAME', 'STATES')
  } else if ( level == '6' ) {
    # 387 features
    #
    # [1] "SHAPE_AREA" "SOURCEFEAT" "AREASQKM"   "METASOURCE" "HUC6"       "SOURCEORIG" "LOADDATE"  
    # [8] "TNMID"      "GEODB_OID"  "OBJECTID"   "AREAACRES"  "GNIS_ID"    "NAME"       "SOURCEDATA"
    #[15] "STATES"     "SHAPE_LENG"
    usefulColumns <- c('AREASQKM', 'HUC6', 'NAME', 'STATES')
  } else if (level == '8' ) {
    # 2300 features
    #
    # [1] "SHAPE_AREA" "SOURCEFEAT" "AREASQKM"   "METASOURCE" "SOURCEORIG" "HUC8"       "LOADDATE"  
    # [8] "TNMID"      "AREAACRES"  "GNIS_ID"    "NAME"       "SOURCEDATA" "STATES"     "SHAPE_LENG"
    usefulColumns <- c('AREASQKM', 'HUC8', 'NAME', 'STATES')
  } else if (level == '10' ) {
    # 18409 features
    # Warning:  Dropping null geometries: 17248 (WBDHU10_01)
    #
    # [1] "SHAPE_AREA" "SOURCEFEAT" "AREASQKM"   "METASOURCE" "SOURCEORIG" "HUMOD"      "LOADDATE"  
    # [8] "TNMID"      "AREAACRES"  "GNIS_ID"    "NAME"       "SOURCEDATA" "STATES"     "HUC10"     
    #[15] "HUTYPE"     "SHAPE_LENG"
    usefulColumns <- c('AREASQKM', 'HUC10', 'NAME', 'STATES')
  } else if (level == '12' ) {
    # 100537 features
    # Warning:  Dropping null geometries: 21407, 21446, 21453, 21869, 21886, 21917, 31625, 31652, 80990, 81132 (WDBHU12_01)
    #
    # [1] "AREAACRES"  "AREASQKM"   "GNIS_ID"    "HUC12"      "HUMOD"      "HUTYPE"     "LOADDATE"  
    # [8] "METASOURCE" "NAME"       "NONCONTR_A" "NONCONTR_K" "OBJECTID"   "SHAPE_AREA" "SHAPE_LEN" 
    #[15] "SOURCEDATA" "SOURCEFEAT" "SOURCEORIG" "STATES"     "TNMID"      "TOHUC"     
    usefulColumns <- c('AREASQKM', 'HUC12', 'NAME', 'STATES')
  } else if (level == '14' ) {
    # NOTE:  On 2015-12-04 it looks like the WBDHU14 file only contains HUCs for southeast Alaska
    # 6865 features
    # Warning: Dropping null geometries: 43, 663, 746, 986, 1289, 1292, 1537, 2189, 2418, 2832, 3720, 3908, 4093, 4248, 4371, 4567, 4652, 4656, 4658, 5528, 5571, 5606, 5691, 6745, 6778 (WBDHU14_01)
    #
    # [1] "AREAACRES"  "AREASQKM"   "GNIS_ID"    "HUC14"      "HUMOD"      "HUTYPE"     "LOADDATE"  
    # [8] "METASOURCE" "NAME"       "NONCONTR_A" "NONCONTR_K" "SHAPE_AREA" "SHAPE_LENG" "SOURCEDATA"
    #[15] "SOURCEFEAT" "SOURCEORIG" "STATES"     "TNMID"     
    usefulColumns <- c('AREASQKM', 'HUC14', 'NAME', 'STATES')
  }
  
  SPDF <- SPDF[,usefulColumns]
  names(SPDF) <- c('area','HUC','HUCName', 'allStateCodes')

  # Change are from km^2 to m^2
  SPDF@data$area <- SPDF@data$area * 1000000
  
  # Group polygons with duplicated hydrologic unit codes
  # NOTE:  The USGS WBD polygons seem to be well organized
  if ( length(SPDF@polygons) != nrow(SPDF@data) ) {
    SPDF <- organizePolygons(SPDF, uniqueID='HUC', sumColumns='area')
  }

  # TODO:  Larger HUCs are centered in the US, while at smaller levels the entire
  # TODO:  HUCs are in foreign countries (ie Canada). Find a way to eliminate smaller
  # TODO:  HUCs whose 'allStateCode' is not a US State

  # Calculate centroids to help add more metadata
  result <- try( {
    centroids <- rgeos::gCentroid(SPDF, byid=TRUE)
    lon <- sp::coordinates(centroids)[,1]
    lat <- sp::coordinates(centroids)[,2]
  }, silent=TRUE)
  
  # NOTE:  This failed for a simplified version of HU10 with:
  # NOTE:
  # NOTE:  Error in createPolygonsComment(p) : 
  # NOTE:    rgeos_PolyCreateComment: orphaned hole, cannot find containing polygon for hole at index 147
  # NOTE:
  # NOTE:  If centroids don't work we'll just default to the center of the bbox for each polygon
  
  if ( class(result)[1] == "try-error" ) {
    cat(paste0('NOTE: rgeos::gCentroid() failed with the following message. Using bbox() to calculate lon and lat.\n'))
    cat(paste0(geterrmessage(),'\n'))
    lon <- rep(as.numeric(NA), nrow(SPDF))
    lat <- rep(as.numeric(NA), nrow(SPDF))
    for (i in 1:nrow(SPDF)) {
      bbox <- bbox(SPDF[i,])
      lon[i] <- mean(bbox[1,])
      lat[i] <- mean(bbox[2,])
    }
  }
  
  # Add more standard columns
  SPDF$longitude <- lon
  SPDF$latitude <- lat  
  SPDF$countryCode <- 'US'
  SPDF$countryName <- 'United States'
  
  #NOTE: this takes quite a long time. 
  suppressWarnings(SPDF$stateCode <- getStateCode(lon, lat, countryCodes=c('US')))
   
  # Hack to change missing stateCodes to the value from allStateCode
  
  for (i in 1:nrow(SPDF)){
    if (is.na(SPDF@data$stateCode[i])){
      SPDF@data$stateCode[i] <- SPDF@data$allStateCode[i]      
    }
    if (stringr::str_length(SPDF@data$stateCode[i]) > 2){
      SPDF@data$stateCode[i] <- substr(SPDF@data$stateCode[i], start=1, stop=2)
    }
  }
  
  
  SPDF$stateName <- codeToState(SPDF$stateCode, SPDF$countryCode)
   
  # Assign a name and save the data
  assign(datasetName,SPDF)
  save(list=c(datasetName),file=paste0(dataDir,"/",datasetName, '.RData'))
  
  return(invisible(datasetName))
}

