### This file is part of 'PGRdup' package for R.

### Copyright (C) 2014, ICAR-NBPGR.
#
# PGRdup is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# PGRdup is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/



#' Double Metaphone Phonetic Algorithm
#' 
#' \code{DoubleMetaphone} converts strings to double metaphone phonetic codes.
#' 
#' An implementation of the Double Metaphone phonetic algorithm in \code{R}. The
#' presence of non-ASCII characters is detected and indicated as a warning.
#' 
#' @param str A character vector whose strings are to be encoded by double 
#'   metaphone algorithm.
#' @return Returns a list with two character vectors of the same length as the
#'   input vector. The first character vector contains the primary double
#'   metaphone encodings, while the second character vector contains the 
#'   alternate encodings.
#' @seealso \code{\link[stringdist]{phonetic}},
#'   \code{\link[RecordLinkage]{phonetics}}
#' @section Acknowledgement: The \code{C} code for the double metaphone
#'   algorithm was adapted from Maurice Aubrey's perl module hosted at the 
#'   \strong{gitpan/Text-DoubleMetaphone} 
#'   \href{https://github.com/gitpan/Text-DoubleMetaphone/blob/master/double_metaphone.c}{public
#'    github library} along with the corresponding 
#'   \href{https://github.com/gitpan/Text-DoubleMetaphone/blob/master/README}{license
#'    information}.
#' @references Philips, L. (2000). The double metaphone search algorithm. C/C++ 
#'   users journal, 18(6), 38-43.
#' @note In case of non-ASCII characters in strings, a warning is issued.
#' @examples
#' # Return the primary and secondary Double Metaphone encodings for a character vector.
#' str1 <- c("Jyothi", "Jyoti")
#' str2 <- c("POLLACHI", "BOLLACHI")
#' DoubleMetaphone(str1)
#' DoubleMetaphone(str2)
#' \dontrun{
#' # Issue a warning in case of non-ASCII characters.
#' str3 <- c("J\xf5geva", "Jogeva")
#' DoubleMetaphone(str3) }
#' @export
#' @useDynLib PGRdup fdouble_metaphone
DoubleMetaphone <- function(str) {
  if (is.character(str) == FALSE) {
    # Check if str is a character vector and stop if not
    stop("str is not a character vector")
  }
  if (any(grepl("NON_ASCII",
                iconv(str, "latin1", "ASCII", sub = "NON_ASCII"))) == TRUE) {
    str <- iconv(str,to = "ASCII//TRANSLIT")
    warning("Non-ASCII characters were encountered.")
  }
  out <- .C("fdouble_metaphone", as.character(str),
            primary = character(length(str)),
            alternate = character(length(str)),
            len = length(str), PACKAGE = "PGRdup")
  out <- out[2:3]
  return(out)
}
