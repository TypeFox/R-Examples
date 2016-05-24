#' Save multiple objects to a file.
#'
#' A simple wrapper for \code{\link[base]{save}}. Understands key = value syntax to save
#' objects using arbitrary variable names. All options of \code{\link[base]{save}},
#' except \code{list} and \code{envir}, are available and passed to
#' \code{\link[base]{save}}.
#'
#' @param file
#'   File to save.
#' @param ... [\code{any}]\cr
#'   Will be converted to an environment and then passed to \code{\link[base]{save}}.
#' @param ascii
#'   See help of \code{\link[base]{save}}.
#' @param version
#'   See help of \code{\link[base]{save}}.
#' @param compress
#'   See help of \code{\link[base]{save}}.
#' @param compression_level
#'   See help of \code{\link[base]{save}}.
#' @param eval.promises
#'   See help of \code{\link[base]{save}}.
#' @param precheck
#'   See help of \code{\link[base]{save}}.
#' @return See help of \code{\link[base]{save}}.
#' @export
#' @examples
#' x = 1
#' save2(y = x, file = tempfile())
save2 = function(file, ..., ascii = FALSE, version = NULL, compress = !ascii,
  compression_level, eval.promises = TRUE, precheck = TRUE) {
  args = tryCatch(as.environment(argsAsNamedList(...)), error = function(e) stopf("Unable to convert to environment (%s)", as.character(e)))
  save(list = ls(args, all.names = TRUE), envir = args, file = file, ascii = ascii, version = version, compress = compress,
       compression_level = compression_level, eval.promises = eval.promises, precheck = precheck)
}
