## ----, eval = FALSE------------------------------------------------------
#  library(MazamaSpatialUtils)
#  
#  # Set data directory to working directory
#  setSpatialDataDir('./SpatialData')
#  initializeSpatialData()

## ----, eval = FALSE------------------------------------------------------
#  convertNaturalEarthAdm1()
#  convertTMWorldBorders()
#  convertTMWorldBordersSimple()
#  convertUSCensusCounties()
#  convertWorldTimezones()

## ----, eval = FALSE------------------------------------------------------
#  loadSpatialData('NaturalEarthAdm1')
#  loadSpatialData('TMWorldBorders')
#  loadSpatialData('TMWorldBordersSimple')
#  loadSpatialData('USCensusCounties')
#  loadSpatialData('WorldTimezones')

## ------------------------------------------------------------------------
library(MazamaSpatialUtils)

lon <- c(-122.3, -73.5, 21.1, 2.5)
lat <- c(47.5, 40.75, 52.1, 48.5)

getCountry(lon, lat)

getCountryCode(lon, lat)

getCountry(lon, lat, allData=TRUE)

## ----, eval=FALSE--------------------------------------------------------
#  # Load states dataset if you haven't already
#  loadSpatialData("NaturalEarthAdm1")
#  
#  # Get which countries the points are in
#  countryCodes <- getCountryCode(lon, lat)
#  
#  # Pass the country codes as an argument to speed everything up
#  getState(lon, lat, countryCodes = countryCodes)
#  
#  getStateCode(lon, lat, countryCodes = countryCodes)
#  
#  # This is a very detailed dataset so we'll grab a few important columns
#  states <- getState(lon, lat, allData=TRUE, countryCodes = countryCodes)
#  states[c('countryCode', 'stateCode', 'stateName')]

## ------------------------------------------------------------------------
# Find the timezones the points are in
getTimezone(lon, lat)

# Find which countries the points are in
countryCodes <- getCountryCode(lon, lat)

# Pass the country codes as an argument to potentially speed things up
getTimezone(lon, lat, countryCodes = countryCodes)

getTimezone(lon, lat, allData=TRUE, countryCodes = countryCodes)

## ----, eval=FALSE--------------------------------------------------------
#  # Load counties dataset if you haven't already
#  loadSpatialData("USCensusCounties")
#  
#  # New dataset of points only in the US
#  stateCodes <- getStateCode(lon,lat)
#  
#  # Optionally pass the state codes as an argument to speed everything up
#  getUSCounty(lon, lat, stateCodes=stateCodes)
#  
#  getUSCounty(lon, lat, allData=TRUE, stateCodes=stateCodes)

## ----, eval=FALSE--------------------------------------------------------
#  library(sp)         # For spatial plotting
#  
#  # Assign timezones polygons an index based on UTC_offset
#  colorIndices <- .bincode(SimpleTimezones@data$UTC_offset, breaks=seq(-12.5,12.5,1))
#  
#  # Color our timezones by UTC_offset
#  plot(SimpleTimezones, col=rainbow(25)[colorIndices])
#  title(line=-3,'Timezone Offsets from UTC')

## ----, eval=FALSE--------------------------------------------------------
#  library(sp)         # For spatial plotting
#  
#  # Read in ISO-encoded oil production and consumption data
#  prod <- read.csv(url('http://mazamascience.com/OilExport/BP_2014_oil_production_bbl.csv'),
#                   skip=6, stringsAsFactors=FALSE, na.strings='na')
#  cons <- read.csv(url('http://mazamascience.com/OilExport/BP_2014_oil_consumption_bbl.csv'),
#                   skip=6, stringsAsFactors=FALSE, na.strings='na')
#  
#  # Only work with ISO-encoded columns of data
#  prodCountryCodes <- names(prod)[ stringr::str_length(names(prod)) == 2 ]
#  consCountryCodes <- names(cons)[ stringr::str_length(names(cons)) == 2 ]
#  
#  # Use the last row (most recent data)
#  lastRow <- nrow(prod)
#  year <- prod$YEAR[lastRow]
#  
#  # Neither dataframe contains all countries so create four categories based on the
#  # amount of information we have:  netExporters, netImporters, exportOnly, importOnly
#  sharedCountryCodes <- intersect(prodCountryCodes,consCountryCodes)
#  net <- prod[lastRow, sharedCountryCodes] - cons[lastRow, sharedCountryCodes]
#  
#  # Find codes associated with each category
#  netExportCodes <- sharedCountryCodes[net > 0]
#  netImportCodes <- sharedCountryCodes[net <= 0]
#  exportOnlyCodes <- setdiff(prodCountryCodes,consCountryCodes)
#  importOnlyCodes <- setdiff(consCountryCodes,prodCountryCodes)
#  
#  # Create a logical 'mask' associated with each category
#  netExportMask <- SimpleCountries@data$countryCode %in% netExportCodes
#  netImportMask <- SimpleCountries@data$countryCode %in% netImportCodes
#  onlyExportMask <- SimpleCountries@data$countryCode %in% exportOnlyCodes
#  onlyImportMask <- SimpleCountries@data$countryCode %in% importOnlyCodes
#  
#  color_export = '#40CC90'
#  color_import = '#EE5555'
#  color_missing = 'gray90'
#  
#  # Base plot (without Antarctica)
#  notAQ <- SimpleCountries@data$countryCode != 'AQ'
#  plot(SimpleCountries[notAQ,],col=color_missing)
#  
#  plot(SimpleCountries[netExportMask,],col=color_export,add=TRUE)
#  plot(SimpleCountries[onlyExportMask,],col=color_export,add=TRUE)
#  plot(SimpleCountries[netImportMask,],col=color_import,add=TRUE)
#  plot(SimpleCountries[onlyImportMask,],col=color_import,add=TRUE)
#  
#  legend('bottomleft',legend=c('Net Exporters','Net Importers'),fill=c(color_export,color_import))
#  title(line=0,paste('World Crude Oil in',year))

