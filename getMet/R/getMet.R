#######################################################################################################
#               Written by Andrew R Sommerlot <andrewrs@vt.edu>, November 2015                        #
#######################################################################################################
#' Gets met data from a specified source and creates model input files in the specified format
#' @param locations - data.frame object or location of csv file with two columns, the first being lattitude and second being longidute in decimal degrees. These locations will be the locations of the collected met data.
#' @param dataSource - Source of met data from predefined source list. Currenly only 'cfsr' is supported.
#' @param outFormat - Format of met data output from predefined source. Currently only 'swat' is supported.
#' @param outDir - Directory where ouput files will be saved
#' @param apiKey - String input of api key for selected data source if required.
#' @return returns specified  met data in specified format
#' @examples
#' \dontrun{
#' locations = data.frame(lat = 38, lon = 79)
#' outDir = "test"
#' getMet(locations=locations, outDir=outDir, dataSource = 'cfsr', outFormat = 'swat')
#' }
#' @export


getMet <-  function(locations, dataSource = 'cfsr', outFormat = 'swat', outDir = getwd(), apiKey = '') {

  if(dataSource != 'cfsr' && dataSource != 'wunderForecast'){
      stop('The getMet function does not yet have support for this data source\n', 'Set dataSource to "cfsr" or "wunderForecast" to continue')
  }
  if(outFormat != 'swat') {
    stop('The getMet function does not yet have support for output formats other than the swat model\n', 'Set outFormat to "swat" to continue')
  }
  if(dataSource == 'cfsr') {
    getSWATcfsr(locations, outDir = outDir)
  }
  if(dataSource == 'wunderForecast' && apiKey == '') {
      stop('Using wunderForecast requires and api key. Get one for free at https://www.wunderground.com/weather/api/d/pricing.html. Select the free version of the ANVIL plan.')
  }
  if(dataSource == 'wunderForecast' && apiKey != '') {
    getSWATwunderForecast(locations, outDir = outDir, apiKey = apiKey)
  }
}


