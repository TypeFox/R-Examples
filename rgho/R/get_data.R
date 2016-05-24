#' Returns GHO Data
#'
#' Given a dimension and a code, returns the corresponding GHO data.
#'
#' Filtering parameters are given as a named list of the form
#' \code{list(COUNTRY = "FRA", ...)}.
#'
#' Query parameters follow the specification described on the WHO website
#' \url{http://apps.who.int/gho/data/node.resources.api}.
#'
#' @param code A GHO code.
#' @param dimension A GHO dimension.
#' @param filter A named list of filtering parameters (see details).
#' @param ... Additional query parameters (see details).
#'
#' @return A \code{data_frame}.
#'
#' @export
#'
#' @examples
#'
#' result <- get_gho_data(
#'   dimension = "GHO",
#'   code = "MDG_0000000001"
#' )
#' print(result, width = Inf)
#'
#'
#' result <- get_gho_data(
#'   dimension = "GHO",
#'   code = "MDG_0000000001",
#'   filter = list(
#'     REGION = "EUR",
#'     YEAR = "2015"
#'   )
#' )
#' print(result, width = Inf)
#'
get_gho_data_ <- function(code, dimension = "GHO", filter = NULL, ...) {

  stopifnot(
    dimension %in% get_gho_dimensions(),
    code %in% get_gho_codes(dimension = dimension)
  )

  get_gho(
    url = build_gho_url(
      dimension = dimension,
      code = code,
      format = "csv",
      filter = filter,
      ...
    )
  ) %>%
    httr::content(type = "text", encoding = "UTF-8") %>%
    readr::read_csv()
}

#' @rdname get_gho_data_
#' @export
get_gho_data <- memoise::memoize(get_gho_data_)
