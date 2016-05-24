#' Taxonomy based cleaning
#'
#' @name taxonomy
#' @param x (data.frame) A data.frame
#' @param name (character) Taxonomic name field Optional. See Details.
#' @param drop (logical) Drop bad data points or not. Either way, we parse
#' out bade data points as an attribute you can access. Default: \code{TRUE}
#' @return Returns a data.frame, with attributes
#' @examples
#' if (requireNamespace("rgbif", quietly = TRUE)) {
#'    library("rgbif")
#'    res <- rgbif::occ_data(limit = 200)$data
#' } else {
#'    res <- sample_data_3
#' }
#'
#' # Remove records where names don't have genus + epithet
#' ## so removes those with only genus and those with no name (NA or NULL)
#' NROW(res)
#' df <- dframe(res) %>% tax_no_epithet(name = "name")
#' NROW(df)
#' attr(df, "name_var")
#' attr(df, "tax_no_epithet")

#' @export
#' @rdname taxonomy
tax_no_epithet <- function(x, name = NULL, drop = TRUE) {
  x <- do_name(x, name)
  noep <- x[!vapply(x$name, function(y) {
    length(strsplit(y, "\\s|_")[[1]])
  }, numeric(1)) >= 2, ]
  if (NROW(noep) == 0) noep <- NA
  if (drop) {
    x <- x[vapply(x$name, function(y) {
      length(strsplit(y, "\\s|_")[[1]])
    }, numeric(1)) >= 2, ]
  }
  row.names(noep) <- NULL
  row.names(x) <- NULL
  structure(x, tax_no_epithet = noep)
}

do_name <- function(x, name) {
  if (is.null(attr(x, "name_var"))) {
    guess_name(x, name)
  } else {
    return(x)
  }
}
