#' Check if the given object is an url string
#'
#' Arguments:
#'  @param s input object to check
#'
#' Returns:
#'  @return TRUE/FALSE indicating whether the input string is a valid URL.
#'
#' @export
#' @references
#' See citation("statfi") 
#' @author Leo Lahti \email{louhos@@googlegroups.com}
#' @examples # is.url("http://aa.px")
#' @keywords utilities
is_url <- function (s) {
  (class(s) == "character" && substr(s,1,7) == "http://")
}

