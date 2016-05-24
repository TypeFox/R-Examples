##' Access to World Bank WDI API
##'
##' A function to extract data from the World Bank API
##'
##' Please refer to \url{http://data.worldbank.org/node/18} for any
##' difference between the country code system. Further details on World
##' Bank classification and methodology are available at
##' \url{http://data.worldbank.org/}.
##'
##' @param indicator The World Bank official indicator name.
##' @param name The new name to be used in the column.
##' @param startDate The start date for the data to begin
##' @param endDate The end date.
##' @param printURL Whether the url link for the data should be printed
##' @param outputFormat The format of the data, can be 'long' or 'wide'.
##' @return A data frame containing the desired World Bank Indicator
##' @export
##'
##' @seealso \code{\link{getFAO}}, \code{\link{getWDItoSYB}},
##' \code{\link{getFAOtoSYB}}
##' @examples
##' ## pop.df = getWDI()
##'

## source:
## http://lamages.blogspot.it/2011/09/setting-initial-view-of-motion-chart-in.html

## TODO (Michael): Investigate why sometimes ISO2 is used and
##                 sometiems ISO3 is used.
getWDI = function(indicator = "SP.POP.TOTL", name = NULL,
                  startDate = 1960, endDate = format(Sys.Date(), "%Y"),
                  printURL = FALSE, outputFormat = "wide"){
    if(is.null(name))
        name = indicator
  url = paste("http://api.worldbank.org/countries/all/indicators/", indicator,
              "?date=", startDate, ":", endDate,
              "&format=json&per_page=30000", sep = "")
  if(printURL)
      print(url)
  wbData = fromJSON(url)[[2]]
  wbData = data.frame(Country = sapply(wbData,
                         function(x) x[["country"]]["value"]),
                       ISO2_WB_CODE= sapply(wbData,
                         function(x) x[["country"]]["id"]),
                       Year = as.integer(sapply(wbData, "[[", "date")),
                       Value = as.numeric(sapply(wbData, function(x)
                         ifelse(is.null(x[["value"]]), NA, x[["value"]]))),
                       stringsAsFactors = FALSE)
  if(outputFormat == "long"){
      wbData$name = name
  } else if(outputFormat == "wide"){
      names(wbData)[4] = name
  }
  wbData
}
