##
#' Generic method to merge two vectors and create a header
#'
#' @param  a Values in the first column.
#' @param  b Values in the second column.
#' @param  header Vector with both column names.
#' @return \code{apa.merge} object; a list consisting of
#' \item{succes}{message in case of an error}
#' \item{data}{the merged data}
#' \item{header}{the header for the merged data}
#' @export
#'
#' @examples
#'
#' # Use apa.merge function
#' apa.merge(a = rnorm(100, mean = 0, sd = 1), b = rnorm(100, mean = 0, sd = 1), header = c("M", "SD"))
##
apa.merge = function(a=NULL, b=NULL, header=NULL) UseMethod("apa.merge")

##
#' Default method to merge two vectors and create a header
#'
#' @param  a Values in the first column.
#' @param  b Values in the second column.
#' @param  header Vector with both column names.
#' @return \code{apa.merge} object; a list consisting of
#' \item{succes}{message in case of an error}
#' \item{data}{the merged data}
#' \item{header}{the header for the merged data}
#' @export
#'
#' @examples
#'
#' # Use apa.merge function
#' apa.merge(a = rnorm(100, mean = 0, sd = 1), b = rnorm(100, mean = 0, sd = 1), header = c("M", "SD"))
##
apa.merge.default = function(a=NULL, b=NULL, header=NULL) {

  est = apaStyleMerge(a, b, header)
  est$call = match.call()
  class(est) = "apa.merge"
  est

}

##
#' Define a print method
#'
#' @param  x A \code{apa.merge} object
#' @export
##
print.apa.merge = function(x, ...) {
  if(x$succes == TRUE) {
    cat("\n")
    cat("Succesfully merged the data.")
    cat("\n\n")
  }
}

# The main function
apaStyleMerge = function(a, b, header) {

  # Initialize function
  options(warn = 0)

  # Check if a valid data frame is supplied
  if ((!is.numeric(a) || !is.numeric(b)) ||
      (is.numeric(a) && length(a) == 0) ||
      (is.numeric(b) && length(b) == 0) ||
      (length(a) != length(b))
    ) {
    error = "Invalid data is supplied."
    warning(error)
    return(list(succes = error))
  }

  # Check if valid headers are supplied
  if((!is.character(header)) || (length(header) != 2)) {
    error = "No valid headers are specified."
    warning(error)
    return(list(succes = error))
  }

  apa.data = paste(sprintf("%3.2f", round(a, digits = 2)), " (", sprintf("%3.2f", round(b, digits=2)), ")", sep = "")
  apa.header = paste(header[1], " (", header[2], ")", sep = "")

  return(list(succes = TRUE, data = apa.data, header = apa.header))

}
