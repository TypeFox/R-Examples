#' Returns GHO Codes for a Given Dimension
#'
#' @param dimension A GHO dimension.
#' @param x For \code{filter}, \code{gho} object to
#'   filter.
#' @param  ... Logical predicates. Multiple conditions are
#'   combined with &.
#'
#' @return GHO codes as a character vector, labels as a
#'   \code{label} attribute and other attributes in a \code{attrs}
#'   attribute, as a \code{data_frame}.
#' @export
#'
#' @examples
#'
#' get_gho_codes(dimension = "GHO")
#'
#' results <- get_gho_codes(dimension = "COUNTRY")
#' filter_attrs(
#'   results,
#'   WHO_REGION_CODE == "EUR"
#' )
#'
get_gho_codes <- function(dimension = "GHO") {
  stopifnot(
    dimension %in% get_gho_dimensions()
  )

  xml_codes <- get_gho(
    url = build_gho_url(dimension = dimension)
  ) %>%
    httr::content() %>%
    xml2::xml_find_all("//Code")

  res <- xml_codes %>%
    xml2::xml_attr("Label")

  labels <- xml_codes[1] %>%
    xml2::xml_find_all("//Code/Display") %>%
    xml2::xml_text()

  build_gho(
    res,
    labels = labels,
    attrs = get_attrs(xml_codes)
  )
}
