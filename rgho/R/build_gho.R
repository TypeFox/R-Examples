#' Check and Build a gho Object
#'
#' @param x A vector of codes.
#' @param labels A vector of labels.
#' @param attrs A \code{data_frame} of attributes.
#'
#' @return A \code{gho} object.
build_gho <- function(x, labels, attrs = NULL) {
  stopifnot(
    is.null(dim(x)),
    is.null(dim(labels)),
    is.character(x),
    is.character(labels),
    "tbl_df" %in% class(attrs) | is.null(attrs),
    length(x) == length(labels),
    if (! is.null(attrs)) nrow(attrs) == length(x) else TRUE
  )
  structure(
    as.vector(x),
    labels = as.vector(labels),
    attrs = attrs,
    class = "gho"
  )
}
