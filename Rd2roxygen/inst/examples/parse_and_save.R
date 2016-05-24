##' Parse the input Rd file and save the roxygen documentation into a file.
##'
##' @param path the path of the Rd file
##' @param file the path to save the roxygen documentation
##' @param usage logical: whether to include the usage section in the output
##' @return a character vector if \code{file} is not specified, or write the vector
##' into a file
##' @export
##' @author Hadley Wickham; modified by Yihui Xie <\url{http://yihui.name}>
parse_and_save <- function(path, file, usage = FALSE) {
  parsed <- parse_file(path)
  output <- create_roxygen(parsed, usage = usage)
  if (missing(file)) output else
	cat(paste(output, collapse = "\n"), file = file)
}
