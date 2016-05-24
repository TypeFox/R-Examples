#' Search the Native Species Resolver
#'
#' @export
#' @param species (character) One or more species names. required.
#' @param country (character) A country name. required.
#' @param stateprovince (character) A state or province name
#' @param countyparish (character) A county or parish name
#' @param ... Further args passed on to \code{\link[httr]{GET}}
#' @references \url{http://bien.nceas.ucsb.edu/bien/tools/nsr/nsr-ws/}
#' @details Currently, only one name is allowed per request. We loop internally
#' over a list of length > 1, but this will still be slow due to only 1
#' name per request.
#'
#' Note that this service can be quite slow.
#' @examples \dontrun{
#' nsr("Pinus ponderosa", "United States")
#' nsr(c("Pinus ponderosa", "Poa annua"), "United States")
#' splist <- c("Pinus ponderosa", "Poa annua", "bromus tectorum", "Ailanthus altissima")
#' nsr(splist, country = "United States")
#'
#' # curl options
#' library("httr")
#' nsr("Pinus ponderosa", "United States", config = verbose())
#' }
nsr <- function(species, country, stateprovince = NULL, countyparish = NULL, ...) {
  tmp <- lapply(species, function(z) {
    args <- orc(list(format = 'json', species = z, country = country,
                     stateprovince = stateprovince, countyparish = countyparish))
    x <- nsr_GET(nsr_base(), args, ...)
    orc(x)
  })
  df <- data.table::rbindlist(tmp, fill = TRUE, use.names = TRUE)
  (df <- data.table::setDF(df))
}

nsr_GET <- function(url, args, ...) {
  x <- httr::GET(url, query = args, ...)
  httr::stop_for_status(x)
  xx <- jsonlite::fromJSON(httr::content(x, "text", encoding = "UTF-8"), FALSE)$nsr_results
  if (length(xx) == 0) NULL else xx[[1]]$nsr_result
}

nsr_base <- function() "http://bien.nceas.ucsb.edu/bien/apps/nsr/nsr_ws.php"
