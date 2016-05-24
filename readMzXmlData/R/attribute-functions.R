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

#' Converts XML attributes.
#'
#' This function converts a XML attribute to \code{character}. 
#'
#' @param attributes \code{character}, XML attributes
#' @param attributeName \code{character}, name XML attributes of attribute to
#'  convert
#' @param required \code{logical}, throw an error if an required attribute is
#'  missing.
#'
#' @return \code{character}, attribute string 
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso \code{\link[readMzXmlData]{.attributeToDouble}},
#'  \code{\link[readMzXmlData]{.attributeTimeToDouble}}
#' @rdname attributeToString
#' @keywords internal
#'
.attributeToString <- function(attributes, attributeName, required=FALSE) {
  a <- unname(attributes[attributeName])
  if (required && (is.null(a) || is.na(a))) {
    stop("Malformed mzXML: attribute ", sQuote(attributeName), " is missing!")
  } else {
    return(a)
  }
}

#' Converts XML attributes.
#'
#' This function converts a XML attribute to \code{double}. 
#'
#' @param attributes \code{character}, XML attributes
#' @param attributeName \code{character}, name XML attributes of attribute to
#'  convert
#' @param required \code{logical}, throw an error if an required attribute is
#'  missing.
#'
#' @return \code{double}
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso \code{\link[readMzXmlData]{.attributeToString}},
#'  \code{\link[readMzXmlData]{.attributeTimeToDouble}}
#' @rdname attributeToDouble
#' @keywords internal
#'
.attributeToDouble <- function(attributes, attributeName, required=FALSE) {
  return(as.double(.attributeToString(attributes=attributes,
                                      attributeName=attributeName,
                                      required=required)))
}

#' Converts XML attributes.
#'
#' This function converts a XML time attribute to \code{double}. 
#'
#' @param attributes \code{character}, XML attributes
#' @param attributeName \code{character}, name XML attributes of attribute to
#'  convert
#' @param required \code{logical}, throw an error if an required attribute is
#'  missing.
#'
#' @return \code{double}
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso \code{\link[readMzXmlData]{.attributeToString}},
#'  \code{\link[readMzXmlData]{.attributeToDouble}}
#' @rdname attributeTimeToDouble
#' @keywords internal
#'
.attributeTimeToDouble <- function(attributes, attributeName, required=FALSE) {
  return(as.double(.grepDouble(.attributeToString(attributes=attributes,
                                                  attributeName=attributeName,
                                                  required=required))))
}

