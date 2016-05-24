#' Returns Available GHO Dimensions
#'
#' @return GHO dimensions as a character vector, and labels as a \code{label} attribute.
#'
#' @export
#'
#' @examples
#'
#' get_gho_dimensions()
#'
get_gho_dimensions <- function() {
  xml_dim <- get_gho(
    url = build_gho_url(dimension = NULL)
  ) %>%
    httr::content() %>%
    xml2::xml_find_all("//Dimension")

  res <- xml_dim %>%
    xml2::xml_attr("Label")

  labels <- xml_dim %>%
    xml2::xml_find_all("//Dimension//Display") %>%
    xml2::xml_text()

  build_gho(
    res, labels = labels
  )
}
