#' @rdname distname-utils
#' 
#' @description
#' \code{check_distname} checks if the distribution specified by
#' the \code{distname} argument is implemented in this package.
#' 
#' @return 
#' \code{check_distname} returns (invisible) that the distribution is implemented, 
#' or throws an error otherwise.
#' @export
#' @examples
#' 
#' check_distname("normal")
#' \dontrun{
#' check_distname("my_great_distribution")
#' }

check_distname <- function(distname) {
  stopifnot(length(distname) == 1,
            is.character(distname))
  if (distname %in% get_distnames()) {
    invisible(paste0("A Lambert W x ", distname, " distribution is implemented."))
  } else if (distname %in% "user") {
    warning("User defined distributions are supported, but no thorough checks are performed. \n ",
            "Please check your results thoroughly.")
  } else {
    stop("Distribution ", distname,  " is not implemented.")
  }
}