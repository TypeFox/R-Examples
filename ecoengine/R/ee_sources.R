#' Ecoengine data sources
#'
#' Returns a full list of data sources supported by the ecoengine
#' @template foptions
#' @export
#' @importFrom httr GET content warn_for_status
#' @importFrom lubridate ymd_hms
#' @return \code{data.frame}
#' @examples 
#' # source_list <- ee_sources()
ee_sources <- function(foptions = list()) {
	# base_url <- "http://ecoengine.berkeley.edu/api/sources/?format=json"
	base_url <- paste0(ee_base_url(), "sources/?format=json")
    data_sources <- GET(base_url, foptions)
    warn_for_status(data_sources)
    ds <- content(data_sources, type = "application/json")
    # sources <- rbind_all(ds$results)
    # Fix to remove ldply
    sources <- ldply(ds$results, function(x) data.frame(x))
    sources$retrieved <- ymd_hms(sources$retrieved)
    sources
}


