#' Writes a Makefile to a file
#'
#' Makefiles, as created by \code{\link{makefile}}, only exist in memory until
#' they are written to a file by this function.
#' @param makefile \code{[MakefileR]}\cr A Makefile, created by
#'   \code{\link{makefile}}
#' @param file_name \code{[character(1)]}\cr Target file name
#' @return The value returned by \code{\link{writeLines}}
#' @export
write_makefile <- function(makefile, file_name) {
  writeLines(format(makefile), file_name)
}
