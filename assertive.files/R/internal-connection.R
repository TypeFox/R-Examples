#' Convert file connections to strings
#' 
#' \code{as.character} method for file connections.
#' @param x A file connection.
#' @param ... Not currently used.
#' @return A string containing the target location of the file connection.
#' @seealso \code{\link[base]{file}}, \code{\link[base]{summary.connection}},
#' \code{\link[base]{as.character}}
#' @examples
#' rprofile <- file.path(R.home("etc"), "Rprofile.site")
#' fcon <- file(rprofile)
#' assertive.base::assert_all_are_true(identical(as.character(fcon), rprofile))
#' close(fcon)
#' @method as.character file
#' @export
as.character.file <- function(x, ...)
{
  # Assertion is to double check that no other package has overwritten the 
  # file class. Want to use assert_is_file_connection(x), but can't have 
  # cyclic dependency.
  if(!inherits(x, c("file", "connection")))
  {
    stop("x is not a file connection.")
  }
  summary(x)$description
}
