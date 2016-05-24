#' Get statistics about BISON downloads.
#'
#' @importFrom dplyr rbind_all
#' @export
#'
#' @param what (character) One of stats (default), search, downnload, or wms. See Details.
#' @param ... Further args passed on to httr::GET. See examples in \code{bison}
#' @return A list of data frame's with names of the list the different data sources
#'
#' @details
#' For the 'what' parameter:
#' \itemize{
#'  \item stats - Retrieve all data provider accumulated statistics.
#'  \item search - Retrieve data provider statistics for BISON searches.
#'  \item download - Retrieve data provider statistics for data downloads from BISON.
#'  \item wms - Retrieve data provider statistics for BISON OGC WMS tile requests.
#' }
#'
#' @examples \dontrun{
#' out <- bison_stats()
#' out <- bison_stats(what='search')
#' out <- bison_stats(what='download')
#' out <- bison_stats(what='wms')
#' out$Arctos
#' out$Harvard_University_Herbaria
#' out$ZooKeys
#' }

bison_stats <- function(what='stats', ...)
{
  what <- match.arg(what, c('stats','search','download','wms'))
  pick <- switch(what, stats='all', search='search', download='download', wms='wms')
  url <- switch(what,
                stats = 'http://bison.usgs.ornl.gov/api/statistics/all',
                search = 'http://bison.usgs.ornl.gov/api/statistics/search',
                download = 'http://bison.usgs.ornl.gov/api/statistics/download',
                wms = 'http://bison.usgs.ornl.gov/api/statistics/wms')

  out <- GET(url, ...)
  stop_for_status(out)
  tt <- content(out, as = "text")
  res <- fromJSON(tt, simplifyVector = FALSE)
  output <- lapply(res$data, function(x){
    df <- rbind_all(lapply(x[[pick]], function(g) data.frame(null2na(g), stringsAsFactors = FALSE)))
    list(name=x$name, resources=do.call(c, x$resources), data=df)
  })
  names(output) <- gsub("\\s", "_", vapply(res$data, "[[", "", "name"))
  return( output )
}

null2na <- function(x){
  x[ sapply(x, is.null, USE.NAMES = FALSE) ] <- NA
  x
}
