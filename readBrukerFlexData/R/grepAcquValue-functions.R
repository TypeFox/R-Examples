## Copyright 2010-2012 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of readBrukerFlexData for R and related languages.
##
## readBrukerFlexData is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## readBrukerFlexData is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with readBrukerFlexData. If not, see <http://www.gnu.org/licenses/>

#' Extracts values from acqu file.
#'
#' This function is a helper function to extract values from acqu file.
#' 
#' @param patternStr \code{character}, pattern to look for
#' @param srcStr \code{character} where to look for \code{patternStr}
#' 
#' @return \code{character} vector of the value given in \code{patternStr}
#'
#' @seealso
#'  \code{\link[readBrukerFlexData]{.grepAcquDoubleValue}},
#'  \code{\link[readBrukerFlexData]{.readAcquFile}}
#' @rdname grepAcquValue
#' @keywords internal
.grepAcquValue <- function(patternStr, srcStr) {
  tmpLine <- grep(pattern=patternStr, x=srcStr, value=TRUE)
  
  ## format e.g.
  ## DATATYPE= CONTINUOUS MASS SPECTRUM
  ## .IONIZATION MODE= LD+
  ## $INSTRUM= <AUTOFLEX>
  
  ## remove front/back pattern environment
  tmpLine <- gsub(pattern="(^.*= *<?)|(>? *$)", replacement="", x=tmpLine)
  
  return(tmpLine)
}

#' Extracts values from acqu file.
#'
#' This function is a helper function to extract double values from acqu file.
#' 
#' @param patternStr \code{character}, pattern to look for
#' @param srcStr \code{character} where to look for \code{patternStr}
#' 
#' @return \code{double} vector of the value given in \code{patternStr}
#'
#' @seealso
#'  \code{\link[readBrukerFlexData]{.grepAcquValue}},
#'  \code{\link[readBrukerFlexData]{.readAcquFile}}
#' @rdname grepAcquDoubleValue
#' @keywords internal
.grepAcquDoubleValue <- function(patternStr, srcStr) {
  strValue <- .grepAcquValue(patternStr, srcStr)
  
  ## replace comma by dot
  strValue <- gsub(",", replacement=".", strValue)
  
  return(as.double(strValue))
}
 
