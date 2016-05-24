#' @encoding UTF-8
#' @title Trim white spaces
#' @description Simply trims spaces from the start, end, and within of a string
#' @param x is a character vector.
#' @param delim is the delimiter, default is white spaces \code{" "}
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' @export
# trim(" Daniel   Marcelino   Silva ")
`trim` <- function(x, delim = " ") {
  gsub("^\\s+|\\s+$", "",
       gsub(sprintf("\\s+[%s]\\s+|\\s+[%s]|[%s]\\s+",
                    delim, delim, delim), delim, x))
}### end -- trim function
NULL
