#' Can some strings be used for column or list element names without problems?
#'
#' @param x [\code{character}]\cr
#'   Character vector to check.
#' @param unique [\code{logical(1)}]\cr
#'   Should the names be unique?
#'   Default is \code{TRUE}.
#' @return [\code{logical}]. One Boolean entry for each string in \code{x}.
#'   If the entries are not unique and \code{unique} is enabled, the first duplicate will
#'   be \code{FALSE}.
#' @export
isValidName = function(x, unique = TRUE) {
  if (!is.character(x))
    x = as.character(x)
  # check that make.names does not change the string (otherwise it would be invalid),
  # names are unique (for e.g. colnames) and stuff like ..1 is disallowed
  x == make.names(x, isTRUE(unique)) & !grepl("^\\.\\.[0-9]$", x)
}
