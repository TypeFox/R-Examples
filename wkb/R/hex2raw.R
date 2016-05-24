# Convert a string hexadecimal representation to a raw vector

#' Convert String Hex Representation to Raw Vector
#'
#' Converts a string hexadecimal representation to a \code{raw} vector.
#'
#' @param hex character string or character vector containing a hexadecimal
#'   representation.
#' @details Non-hexadecimal characters are removed.
#' @return A \code{\link[base]{raw}} vector.
#'
#'   The return value is a \code{list} of \code{raw} vectors when the argument
#'   \code{hex} contains more than one hexadecimal representation.
#' @examples
#' # create a character string containing a hexadecimal representation
#' hex <- "0101000000000000000000f03f0000000000000840"
#'
#' # convert to raw vector
#' wkb <- hex2raw(hex)
#'
#'
#' # create a character vector containing a hexadecimal representation
#' hex <- c("01", "01", "00", "00", "00", "00", "00", "00", "00", "00", "00",
#'          "f0", "3f", "00", "00", "00", "00", "00", "00", "08", "40")
#'
#' # convert to raw vector
#' wkb <- hex2raw(hex)
#'
#'
#' # create vector of two character strings each containing a hex representation
#' hex <- c("0101000000000000000000f03f0000000000000840",
#'          "010100000000000000000000400000000000000040")
#'
#' # convert to list of two raw vectors
#' wkb <- hex2raw(hex)
#' @seealso \code{raw2hex} in package
#'   \href{http://cran.r-project.org/package=PKI}{\pkg{PKI}}, \code{\link{readWKB}}
#' @export
hex2raw <- function(hex) {
  if(!(is.character(hex) || (is.list(hex) &&
    all(vapply(X = hex, FUN = is.character, FUN.VALUE = logical(1)))))) {
    stop("hex must be a character string or character vector")
  }
  if(is.list(hex) || (length(hex) > 1 &&
     all(vapply(X = hex, FUN = nchar, FUN.VALUE = integer(1)) > 2))) {
    lapply(hex, .hex2raw)
  } else {
    .hex2raw(hex)
  }
}

.hex2raw <- function(hex) {
  hex <- gsub("[^0-9a-fA-F]", "", hex)
  if(length(hex) == 1) {
    if(nchar(hex) < 2 || nchar(hex) %% 2 != 0) {
      stop("hex is not a valid hexadecimal representation")
    }
    hex <- strsplit(hex, character(0))[[1]]
    hex <- paste(hex[c(TRUE, FALSE)], hex[c(FALSE, TRUE)], sep = "")
  }
  if(!all(vapply(X = hex, FUN = nchar, FUN.VALUE = integer(1)) == 2)) {
    stop("hex is not a valid hexadecimal representation")
  }
  as.raw(as.hexmode(hex))
}
