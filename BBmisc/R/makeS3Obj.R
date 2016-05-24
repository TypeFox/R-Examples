#' Simple constructor for S3 objects based on lists.
#'
#' Simple wrapper for \code{as.list} and \code{\link{setClasses}}.
#'
#' @param classes [\code{character}]\cr
#'   Class(es) for constructed object.
#' @param ... [any]\cr
#'   Key-value pairs for class members.
#' @return Object.
#' @export
#' @examples
#' makeS3Obj("car", speed = 100, color = "red")
makeS3Obj = function(classes, ...) {
  assertCharacter(classes, min.len = 1L, any.missing = FALSE)
  setClasses(list(...), classes = classes)
}
