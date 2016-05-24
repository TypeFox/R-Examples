## Copyright 2012 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of MALDIquantForeign for R and related languages.
##
## MALDIquantForeign is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## MALDIquantForeign is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with MALDIquantForeign. If not, see <http://www.gnu.org/licenses/>

#' This method creates a filename for a \code{\linkS4class{AbstractMassObject}}
#' object.
#'
#' @param x a \code{\linkS4class{AbstractMassObject}} object
#' @param fileExtension file type (e.g. "txt", "pdf", ...)
#' @return filename
#'
#' @seealso \code{\linkS4class{AbstractMassObject}}
#' @aliases .composeFilename,AbstractMassObject-method
#' @aliases .composeFilename,list-method
#' @docType methods
#' @keywords internal
#' @noRd
setMethod(f=".composeFilename",
  signature=signature(x="AbstractMassObject"),
  definition=function(x, fileExtension="csv") {

  if (!is.null(metaData(x)$fullName)) {
    if (length(metaData(x)$fullName) > 1) {
      filename <- paste0(metaData(x)$fullName, collapse="_")
    } else {
      filename <- metaData(x)$fullName
    }
  } else {
    filename <- .withoutFileExtension(metaData(x)$file[1])
  }

  filename <- paste(filename, fileExtension, sep=".")

  return(filename)
})

setMethod(f=".composeFilename",
  signature=signature(x="list"),
  definition=function(x, fileExtension="csv") {

  stopifnot(MALDIquant:::.isMassObjectList(x))

  filenames <- unlist(lapply(x, .composeFilename, fileExtension=fileExtension))

  .uniqueBaseFilenames(filenames, fileExtension)
})
