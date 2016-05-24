#' Download full list of ISO-3166 alpha 3 country code.
#'
#' Function returns a data frame with country name and  \code{iso3c} code which is required by the \code{\link{landings}} function to return country specific data
#'
#' @param  foptions additional curl options
#' @export
#' @return data.frame
#' @importFrom httr GET content stop_for_status
#' @importFrom data.table rbindlist
#' @examples \dontrun{
#' of_country_codes()
#'}
of_country_codes <- function(foptions = list()) {
    url <- "http://openfisheries.org/api/landings/countries.json"
    countries_call <- GET(url, foptions)
    stop_for_status(countries_call)
    countries<- content(countries_call)
    countries <- data.frame(rbindlist(countries), stringsAsFactors = FALSE)
    return(countries)
}

#' country_codes
#' 
#' Function has been deprecated. Now replaced by \code{of_country_codes}
#' @export
#' @rdname country_codes-deprecated
country_codes <- function() {
  .Deprecated(new = "of_country_codes", package = "rfisheries", msg = "This function is deprecated, and will be removed in a future version. See ?of_country_codes")
}