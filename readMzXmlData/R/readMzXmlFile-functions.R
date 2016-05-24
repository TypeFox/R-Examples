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

#' Reads mass spectrometry data in mzXML format.
#'
#' Reads mass spectrometry data in mzXML format defined in
#' \url{http://tools.proteomecenter.org/mzXMLschema.php}.
#'
#' @param mzXmlFile \code{character}, path to \emph{mzXML} file which should
#'  be read.
#' @param removeMetaData \code{logical}, to save memory metadata could be
#'  deleted.
#' @param verbose \code{logical}, verbose output?
#'
#' @return A list of spectra and metadata.
#' \itemize{
#'  \item{\code{spectrum$mass}: }{A vector of calculated mass.}
#'  \item{\code{spectrum$intensity}: }{A vector of intensity values.}
#'  \item{\code{metaData}: }{A list of metaData depending on read spectrum.}
#' }
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso \code{\link{readMzXmlDir}},
#' \code{\link[MALDIquantForeign]{importMzXml}}
#' @references Definition of \code{mzXML} format:
#' \url{http://tools.proteomecenter.org/mzXMLschema.php}
#'
#' @keywords IO
#' @rdname readMzXmlFile
#' @export
#' @examples
#'
#' ## load library
#' library("readMzXmlData")
#'
#' ## get examples directory
#' exampleDirectory <- system.file("Examples", package="readMzXmlData")
#'
#' ## read example spectrum
#' spec <- readMzXmlFile(file.path(exampleDirectory, "A1-0_A1.mzXML"))
#'
#' ## print metaData
#' print(spec$metaData)
#'
#' ## plot spectrum
#' plot(spec$spectrum$mass, spec$spectrum$intensity, type="l")
#'
readMzXmlFile <- function(mzXmlFile, removeMetaData=FALSE, verbose=FALSE) {
  scans <- .readMzXmlFile(mzXmlFile=mzXmlFile, removeMetaData=removeMetaData,
                          verbose=verbose);

  if (length(scans) <= 1) {
    return(scans[[1]]);
  } else {
    return(scans);
  }
}

#' Reads mass spectrometry data in mzXML format.
#'
#' Internal function.
#
#' @param mzXmlFile \code{character}, path to \emph{mzXML} file which should
#'  be read.
#' @param removeMetaData \code{logical}, to save memory metadata could be
#'  deleted.
#' @param verbose \code{logical}, verbose output?
#'
#' @return A list of spectra and metadata.
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso \code{\link[readMzXmlData]{readMzXmlFile}}
#' @rdname readMzXmlFile-internal
#' @keywords internal
#'
.readMzXmlFile <- function(mzXmlFile, removeMetaData=FALSE, verbose=FALSE) {
  if (verbose) {
    message("Reading spectrum from ", sQuote(mzXmlFile), " ...")
  }

  if (!file.exists(mzXmlFile)) {
    stop("File ", sQuote(mzXmlFile), " doesn't exists!")
  }

  if (file.info(mzXmlFile)$isdir) {
    stop("Not a mzXML file! ", sQuote(mzXmlFile), " is a directory.")
  }

  ## try to get absolute file path
  mzXmlFile <- normalizePath(mzXmlFile)

  ## read file
  s <- .parseMzXml(file=mzXmlFile, verbose=verbose)

  spectra <- lapply(s$scan, function(x, globalS=s) {
    scan <- list()
    scan$spectrum <- x$spectrum
    scan$metaData$file <- mzXmlFile

    if (!removeMetaData) {
      scan$metaData <- c(scan$metaData, globalS$metaData, x$metaData)
    }
    return(scan)
  })

  return(spectra)
}

