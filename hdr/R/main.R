create_url <- function(params) {

  base_url <- "http://ec2-52-1-168-42.compute-1.amazonaws.com/version/1/"

  params <- Filter(Negate(is.null), params)
  params <- lapply(params, paste0, collapse = ",")

  query <- paste(names(params), unlist(params), sep = "/", collapse = "/")

  paste0(base_url, query)
}

make_request <- function(url) {

  response <- httr::GET(url)
  httr::stop_for_status(response)
  cont <- httr::content(response)

  dfs <- lapply(cont$indicator_value, function(x) {
    structure(data.frame(x, stringsAsFactors = FALSE),
              names = c("iso3c", "id", "year", "value"))
  })

  df <- do.call(rbind, dfs)

  df[["year"]] <- as.integer(df[["year"]])
  df[["value"]] <- as.numeric(df[["value"]])

  df_ind <- data.frame(id = names(cont$indicator_name),
                       indicator = unlist(cont$indicator_name),
                       stringsAsFactors = FALSE)

  df_country <- data.frame(iso3c = names(cont$country_name),
                           country = unlist(cont$country_name),
                           stringsAsFactors = FALSE)

  df <- merge(df, df_ind)
  df <- merge(df, df_country)

  class(df) <- c("tbl_df", "tbl", "data.frame")

  df[, c("id", "indicator", "iso3c", "country", "year", "value")]
}

#' Fetch data from the UNDP Human Development Report
#'
#' @param indicator Numerical or character vector with the indicator id (see details)
#' @param country Character vector
#' @param year Numerical vector (see details for which years are available)
#'
#' The function fetches data from the \href{http://hdr.undp.org}{United Nations
#' Development Programme Human Development Report} API.
#'
#' A dimension can be left as NULL (the default) to get all data for that
#' dimension. The package includes a data frame (\code{hdr_indicators}) with the
#' IDs and human-readable names of the indicators.
#'
#' If the year parameter is not left as NULL, it must be on of the following:
#' 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2011, 2012, 2013.
#'
#' @return A data frame`F
#' @export
#'
#' @examples
#' # Get the Human Development Index for Germany in 2013
#' df <- get_data(indicator = 137506, country = "DEU", year = 2013)
#' head(df)
#'
#' # Leave a dimension as NULL (default) to get all values for that dimension
#' # e.g. all countries and all year for a specific indicator:
#' df <- get_data(103606)
#' head(df)
#'
#' # A data frame with id and indicator names
#' head(hdr_indicators)
get_data <- function(indicator = NULL, country = NULL, year = NULL) {

  if(!is.null(year) && !year %in% c(seq(1980, 2010, 5), 2011:2013)) {
    stop("Year must be one of: 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2011, 2012, 2013")
  }

  params <- list(indicator_id = indicator, country_code = country, year = year)

  url <- create_url(params)

  make_request(url)
}
