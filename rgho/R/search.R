#' Search Codes or Dimensions Labels
#'
#' @param x A character string to find in GHO labels.
#' @param dimension A GHO dimension where codes will be searched.
#' @param gho An object from \code{\link{get_gho_dimensions}} or
#'  \code{\link{get_gho_codes}}.
#'
#' @return A vector of dimensions or codes whose labels mach the search string.
#'  Labels are given as as a \code{label} attribute.
#'
#' @export
#'
#' @examples
#'
#' search_codes("neonatal", dimension = "GHO")
#'
#' result <- get_gho_codes(dimension = "REGION")
#' search_gho(result, "asia")
#'
search_gho <- function(gho, x) {
  pos <- grep(
    tolower(x),
    tolower(attr(gho, "labels")),
    useBytes = TRUE,
    fixed = TRUE
  )

  if (! is.null(attr(gho, "attrs"))) {
    attrs <- attr(gho, "attrs")[pos, ]
  } else {
    attrs <- NULL
  }

  build_gho(
    gho[pos],
    labels = attr(gho, "labels")[pos],
    attrs = attrs
  )
}

#' @rdname search_gho
#' @export
search_dimensions <- function(x) {
  search_gho(get_gho_dimensions(), x)
}

#' @rdname search_gho
#' @export
search_codes <- function(x, dimension = "GHO") {
  search_gho(get_gho_codes(dimension), x)
}
