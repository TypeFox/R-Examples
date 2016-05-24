#' @docType package
#' @name MazamaSpatialUtils
#' @title Mazama Science spatial data and utility functions.
#' @description This package contains code to convert various spatial datasets into .RData files 
#' with uniformly named identifiers including:
#' \itemize{
#' \item{ countryCode -- ISO 3166-1 alpha-2}
#' \item{ countryName -- Country name}
#' \item{ stateCode -- ISO 3166-2 alpha-2}
#' \item{ timezone -- Olson timezone}
#' \item{ longitude -- degrees East}
#' \item{ latitude -- degrees North}
#' \item{ area -- m^2}
#' }
#' The parameters listed above will be found in the @@data slot of each spatial dataset whose source
#' data has an equivalent field. The only field guaranteed to exist in every dataset is \code{countryCode}.
#' 
#' The following additional standards are applied during the data conversion process:
#' \itemize{
#' \item{ all spatial data are converted to a purely geographic projection (\code{CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"}) }
#' \item{ no duplicated rows in the dataframe (conversion to \strong{multi-}polygons) }
#' \item{ lowerCamelCase, human readable names replace original parameter names }
#' \item{ redundant, software-internal or otherwise unuseful data columns may be dropped }
#' \item{ parameters may be added to the @@data dataframe }
#' \item{ latitude and longitude of polygon centroids may be added }
#' }
#' 
#' Utility functions allow users to determine the country, state, county and timezones
#' associated with a set of locations, e.g. environmental monitoring sites.
#' 
#' The uniformity of identifiers in the spatial datasets also makes it easy to generate maps
#' with data from any dataset that uses standard ISO codes for countries or states.
#' 
#' \strong{History}
#' 
#' version 0.4.2 -- patch
#' \itemize{
#'   \item{Added \code{encoding} argument to converLayer().}
#'   \item{Modified convertUSCensusCounties() to use \code{encoding='latin1'}.}
#' }
#'
#' version 0.4.1 -- patch
#' \itemize{
#'   \item{Added \code{useBuffering} to get-Sate,CountryTimezone functions.}
#' }
#'
#' version 0.3.2 -- patch
#' \itemize{
#'   \item{getSpatialData() no longer fails on invliad/missing locations, now returns dataframe rows with all NA.}
#' }
#'
#'#' version 0.3.1 -- addition of buffered search and WorldEEZ polygons
#' \itemize{
#'   \item{Updated included datasets to use \code{"+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"}.}
#'   \item{Addition of buffered search so that locations can find nearby polygons.}
#'   \item{Addition of convertWorldEEZ() function.}
#' }
#' 
#' version 0.2.4 -- patch
#' \itemize{
#'   \item{Updated default projection from \code{"+proj=longlat"} to \code{"+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"}
#'   to support libproj >= 4.9.1}
#' }
#' 
#' version 0.2.3 -- patch
#' \itemize{
#'   \item{Removed unneeded test that fails with sp version 1.1-0.}
#' }
#' 
#' version 0.2.2 -- minor tweaks to 0.2.1
#' \itemize{
#'   \item{User specification of \code{SpatialDataDir} is now required.}
#'   \item{Minor documentation improvements.}
#' }
#' 
#' version 0.2.1 -- addition of GADM and USGS HUC8
#' \itemize{
#'   \item{Converts for GADM administrative boundaries and and USGS watershed datasets.}
#'   \item{Addition of code-name, name-code and code-code conversion utilities.}
#'   \item{Addition of organizePolygons() function.}
#' }
#' 
#' version 0.1 -- initial release
NULL


#' @docType data
#' @keywords datasets
#' @name SimpleCountries
#' @title World Country Polygons
#' @format A SpatialPolygonsDataFrame with 246 elements and 7 columns of data.
#' @description SimpleCountries is a simplified world borders dataset suitable for global maps
#' and quick spatial searches. This dataset is distributed with the package and is used by
#' default whenever a dataset with country polygons is required.
#' @details This dataset is equivalent to TMWorldBordersSimple but with fewer columns of data.
#' @seealso convertTMWorldBordersSimple
NULL


#' @docType data
#' @keywords datasets
#' @name SimpleTimezones
#' @title World Timezone Polygons
#' @format A SpatialPolygonsDataFrame with 1106 elements and 6 columns of data.
#' @description SimpleTimezones is a simplified world timezones dataset suitable for global maps
#' and quick spatial searches. This dataset is distributed with the package and is used by
#' default whenever a dataset with timezone polygons is required.
#' @details This dataset is a simplified version of WorldTimezones. It was simplified with
#' \url{http://mapshaper.org}.
#' @seealso convertWorldTimezones
NULL


# ----- Internal Package State -------------------------------------------------

spatialEnv <- new.env(parent = emptyenv())
spatialEnv$dataDir <- NULL

#' @docType data
#' @keywords environment
#' @name SpatialDataDir
#' @title Directory for Spatial Data
#' @format Absolute path string.
#' @description This package maintains an internal directory location which users can set
#' using \code{setSpatialDataDir()}. All package functions use this directory whenever datasets
#' are created or loaded.
#' 
#' The default setting when the package is loaded is \code{getwd()}.
#' @seealso getSpatialDataDir
#' @seealso setSpatialDataDir
NULL

#' @keywords environment
#' @export
#' @title Get Package Data Directory
#' @description Returns the package data directory where spatial data is located.
#' @return Absolute path string.
#' @seealso dataDir
#' @seealso setSpatialDataDir
getSpatialDataDir <- function() {
  if (is.null(spatialEnv$dataDir)) {
    stop('No data directory found. Please set a data directory with setSpatialDataDir("YOUR_DATA_DIR").',call.=FALSE)
  } else {
    return(spatialEnv$dataDir)    
  }
}

#' @keywords environment
#' @export
#' @title Set Package Data Directory
#' @param dataDir directory where spatial datasets are created
#' @description Sets the package data directory where spatial data is located.
#' If the directory does not exist, it will be created.
#' @return Silently returns previous value of data directory.
#' @seealso SpatialDataDir
#' @seealso getSpatialDataDir
setSpatialDataDir <- function(dataDir) {
  old <- spatialEnv$dataDir
  dataDir <- path.expand(dataDir)
  tryCatch({
    if (!file.exists(dataDir)) dir.create(dataDir)
    spatialEnv$dataDir <- dataDir
  }, warning = function(warn) {
    warning("Invalid path name.")
  }, error   = function(err) {
    stop(paste0("Error in setSpatialDataDir(",dataDir,")."))
  })     
  return(invisible(old))
}


# ----- Code <-> Name conversion functions  ------------------------------------

#' @keywords conversion
#' @export
#' @title Convert Between ISO2 and ISO3 Country Codes
#' @param countryCodes vector of country codes to be converted
#' @description Converts a vector of ISO 3166-1 alpha-2 codes to the corresponding ISO 3166-1 alpha-3 codes or vice versa.
#' @return A vector of ISO country codes
codeToCode <- function(countryCodes) {
  countryTable <- MazamaSpatialUtils::SimpleCountries@data
  nonMissingCountryCodes <- countryCodes[!is.na(countryCodes)]
  if ( all(stringr::str_length(nonMissingCountryCodes) == 2) ) {
    # Create a vector of ISO3 identified by countryCode
    allISO3 <- countryTable$ISO3
    names(allISO3) <- countryTable$countryCode
    return(as.character(allISO3[countryCodes]))
  } else if ( all(stringr::str_length(nonMissingCountryCodes) == 3) ) {
    # Create a vector of ISO2 identified by ISO3
    allISO2 <- countryTable$countryCode
    names(allISO2) <- countryTable$ISO3
    return(as.character(allISO2[countryCodes]))    
  } else {
    stop('countryCodes must be either all ISO 3166-1 alpha-2 or all ISO 3166-1 alpha-3', call.=FALSE)
  }
}

#' @keywords conversion
#' @export
#' @title Convert Country Codes to Country Names
#' @param countryCodes vector of country codes to be converted
#' @description Converts a vector of ISO 3166-1 alpha-2 codes to the corresponding English names.
#' @return A vector of English country names or NA.
codeToCountry <- function(countryCodes) {
  countryTable <- MazamaSpatialUtils::SimpleCountries@data
  # Create a vector of countryNames identified by countryCode
  allNames <- countryTable$countryName
  names(allNames) <- countryTable$countryCode
  return(as.character(allNames[countryCodes]))
}

#' @keywords conversion
#' @export
#' @title Convert Country Names to Country Codes
#' @param countryNames vector of country names to be converted
#' @description Converts a vector of English country names to the corresponding ISO 3166-1 alpha-2 codes.
#' @return A vector of ISO 3166-1 alpha-2 codes or NA.
countryToCode <- function(countryNames) {
  countryTable <- MazamaSpatialUtils::SimpleCountries@data
  # Create a vector of countryCodes identified by countryName
  allCodes <- countryTable$countryCode
  names(allCodes) <- countryTable$countryName
  return(as.character(allCodes[countryNames]))
}

#' @keywords conversion
#' @export
#' @title Convert State Codes to State Names
#' @param stateCodes vector of state codes to be converted
#' @param countryCodes ISO-3166-1 alpha-2 country codes the state might be found in
#' @param dataset name of dataset containing state-level identifiers
#' @description Converts a vector of ISO 3166-2 alpha-2 state codes to the corresponding English names.
#' @details For this function to work, you must first run \code{initializeSpatialData()} to
#' download, convert and install the necessary spatial data.
#' @return A vector of English state names or NA.
#' @seealso convertNaturalEarthAdm1
codeToState <- function(stateCodes, countryCodes=NULL,
                        dataset='NaturalEarthAdm1') {
  
  # Sanity check
  if (!exists(dataset)) {
    stop('Missing database. Please loadSpatialData("',dataset,'")',call.=FALSE)
  }
  SPDF <- get(dataset)
  # Remove NA state codes
  stateTable <- SPDF@data[!is.na(SPDF@data$stateCode),]
  # Filter by countryCodes to make searching faster
  if (!is.null(countryCodes)) stateTable <- stateTable[stateTable$countryCode %in% countryCodes,]
  # Create a vector of state names identified by state code
  allStates <- stateTable$stateName
  names(allStates) <- stateTable$stateCode
  return(as.character(allStates[stateCodes]))
}

#' @keywords conversion
#' @export
#' @title Convert State Names to State Codes
#' @param stateNames state names to be converted
#' @param countryCodes ISO 3166-2 alpha-2 country codes the state might be found in
#' @param dataset name of dataset containing state-level identifiers
#' @description Converts a vector of state names to an ISO 3166-2 two character state codes.
#' @details For this function to work, you must first run \code{initializeSpatialData()} to
#' download, convert and install the necessary spatial data.
#' @return A vector of ISO 3166-2 codes or NA.
#' @seealso convertNaturalEarthAdm1
stateToCode <- function(stateNames, countryCodes=NULL,
                        dataset='NaturalEarthAdm1') {
  # Sanity check
  if (!exists(dataset)) {
    stop('Missing database. Please loadSpatialData("',dataset,'")',call.=FALSE)
  }
  SPDF <- get(dataset)
  # Remove NA state codes
  stateTable <- SPDF@data[!is.na(SPDF@data$stateCode),]
  # Filter by countryCodes to make searching faster
  if (!is.null(countryCodes)) stateTable <- stateTable[stateTable$countryCode %in% countryCodes,]
  # Create a vector of state codes identified by name
  allCodes <- stateTable$stateCode
  names(allCodes) <- stateTable$stateName
  return(as.character(allCodes[stateNames]))
}


# ----- Initialization----------------------------------------------------------

#' @export
#' @title Install Core Datasets
#' @description Four core datasets can be installed to enhance the base the functionality
#' in \pkg{MazamaSpatialUtils}. Running \code{initializeSpatialData()} will
#' install these datasets in the the directory specified by \code{setSpatialDataDir()}.
#' 
#' The core datastes are:
#' \itemize{
#' \item{\code{TMWorldBorders} -- high resolution country polygons (higher resolution than \code{SimpleCountries})}
#' \item{\code{NaturalEarthAdm1} -- state/province polygons throughout the world}
#' \item{\code{USCensusCounties} -- county polygons in the United States}
#' \item{\code{WorldTimezones} -- high resolution timezone polygons (higher resolution than \code{SimpleTimezones})}
#' }
#' @return Nothing.
#' @seealso setSpatialDataDir
initializeSpatialData <- function() {
  # Install high resolution and load Country, State and Timezone datasets
  installSpatialData('TMWorldBorders')
  loadSpatialData('TMWorldBorders')
  installSpatialData('NaturalEarthAdm1')
  loadSpatialData('NaturalEarthAdm1')
  installSpatialData('USCensusCounties')
  loadSpatialData('USCensusCounties')
  installSpatialData('WorldTimezones')
  loadSpatialData('WorldTimezones')
}

