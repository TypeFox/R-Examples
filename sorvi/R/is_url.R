# To call in the statistician after the experiment is done may be no more
# than asking him to perform a post-mortem examination: he may be able to
# say what the experiment died of. ~ Sir Ronald Aylmer Fisher

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
#' See citation("sorvi") 
#' @author Leo Lahti \email{louhos@@googlegroups.com}
#' @examples is_url("http://aa.px")
#' @keywords utilities
is_url <- function (s) {
  (class(s) == "character" && substr(s,1,7) == "http://")
}



