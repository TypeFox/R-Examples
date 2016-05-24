## Copyright 2011-2012 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of readMzXmlData for R and related languages.
##
## readMzXmlData is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## readMzXmlData is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with readMzXmlData. If not, see <http://www.gnu.org/licenses/>

#' Pattern matching.
#'
#' This function returns a substring selected by a regular expression.
#'
#' @param pattern \code{character}, regexp pattern
#' @param x \code{character}, text
#'
#' @return \code{character}, matched substring
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso \code{\link[readMzXmlData]{.grepNumber}}
#' @rdname grepSubString
#' @keywords internal
#'
.grepSubString <- function(pattern, x) {
  rx <- regexpr(pattern=pattern, text=x)
  str <- substr(x=x, start=rx, stop=(rx + (attr(rx, "match.length") - 1) ))
  return(str)
}

#' Pattern matching.
#'
#' This function returns \code{numeric} values found in a \code{character}.
#'
#' @param x \code{character}, text
#'
#' @return \code{numeric}, matched value
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso \code{\link[readMzXmlData]{.grepDouble}}
#'  \code{\link[readMzXmlData]{.grepSubString}}
#' @rdname grepNumeric
#' @keywords internal
#'
.grepNumber <- function(x) {
  return(.grepSubString(pattern="[0-9]+\\.?[0-9]*", x=x))
}

#' Pattern matching.
#'
#' This function returns \code{double} values found in a \code{character}.
#'
#' @param x \code{character}, text
#'
#' @return \code{double}, matched value
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso \code{\link[readMzXmlData]{.grepNumber}}
#'  \code{\link[readMzXmlData]{.grepSubString}}
#' @rdname grepNumeric
#' @keywords internal
#'
.grepDouble <- function(x) {
  return(as.double(.grepNumber(x)))
}

