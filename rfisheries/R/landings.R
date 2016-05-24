#' Returns landings data from the openfisheries API
#'
#' The function returns aggregate landings data if no parameters are supplied. One could get country or species-specific data by specifying either one of those options. Country must be provided as the \code{iso3c} code and species must be supplied as a3_code. Supporting functions \code{country_codes} and \code{species_codes} provide that data and can be combined to return data for multiple countries or species.
#' @param country Default is \code{NA}. Download country specific data by specifying the ISO-3166 alpha 3 country code.
#' @param  species Default is \code{NA}. Download species specific data by specifying the three-letter ASFIS species code
#' @param foptions additional optional parameters
#' @export
#' @import rjson
#' @importFrom httr GET content stop_for_status
#' @importFrom data.table rbindlist
#' @return data.frame
#' @examples \dontrun{
#' of_landings()
#' # Landings by country
#' of_landings(country = 'CAN')
#' #landings by species
#' of_landings(species = 'COD')
#'}
of_landings <- function(country = NA, species = NA, foptions = list()) {
    if (!is.na(country) && !is.na(species))
        stop("Specify country or species but not both", call. = FALSE)
    if (is.na(country) && is.na(species)) {
        url <- "http://openfisheries.org/api/landings.json"
    } else if (!is.na(country) && is.na(species)) {
        url <- paste0("http://openfisheries.org/api/landings/countries/",
            country, ".json")
    } else {
        url <- paste0("http://openfisheries.org/api/landings/species/",
            species, ".json")
    }
    landings_call <- GET(url, foptions)
    stop_for_status(landings_call)
    landings_data_JSON <- content(landings_call)
    if (length(landings_data_JSON) == 0) {
        landings_data <- data.frame()
    } else {
    landings_data <- data.frame(rbindlist(landings_data_JSON))
    # Add the species as a column to avoid ambguity
    if(!is.na(species))  {
        landings_data <- cbind(landings_data, species)
        landings_data$species <- as.character(landings_data$species)
      }
      # Do the same with the country.
      if (!is.na(country)) {
        landings_data <- cbind(landings_data, country)
        landings_data$country <- as.character(landings_data$country)
      }
    }

    if (nrow(landings_data) == 0) {
      # stop("No data found", call. = FALSE)
      NULL
    } else {
      landings_data
    }
}


#' @rdname country_codes-deprecated
#' @export
landings <- function()
{
  .Deprecated(new = "of_landings", package = "rfisheries", msg = "This function is deprecated, and will be removed in a future version. See ?of_landings")
}

