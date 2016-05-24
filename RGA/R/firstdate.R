#' @title First Date
#'
#' @description Get the first date with available data.
#'
#' @param profileId character. Google Analytics profile ID. Can be obtained using the \code{\link{list_profiles}} or via the web interface Google Analytics.
#' @param token \code{\link[httr]{Token2.0}} class object with a valid authorization data.
#'
#' @return Start date of collecting the Google Analytics statistics.
#'
#' @family Reporting API
#'
#' @examples
#' \dontrun{
#' authorize()
#' first_date <- firstdate(profileId = "profile_id")
#' ga_data <- get_ga("profile_id", start.date = first_date, end.date = "today",
#'                   metrics = "ga:sessions", dimensions = "ga:source,ga:medium",
#'                   sort = "-ga:sessions")
#' }
#'
#' @include ga.R
#' @export
firstdate <- function(profileId, token) {
    res <- suppressWarnings(
        get_ga(profileId = profileId, start.date = "2005-01-01", end.date = "today",
               metrics = "ga:sessions", dimensions = "ga:date", filters = "ga:sessions>0",
               max.results = 1L, token = token)
    )
    return(res$date)
}
