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

#' Reads recursively mass spectrometry data in mzXML format.
#'
#' Reads recursively all mass spectrometry data in mzXML format in a specified
#' directory.
#'
#' @details See \code{\link{readMzXmlFile}}.
#'
#' @param mzXmlDir \code{character}, path to \emph{directory} which should
#'  be read recursively.
#' @param removeCalibrationScans \code{logical}, if \code{TRUE} all scans in
#'  directories called \dQuote{[Cc]alibration} will be ignored.
#' @param removeMetaData \code{logical}, to save memory metadata could be
#'  deleted.
#' @param rewriteNames \code{logical}, if \code{TRUE} all list elements get
#'  an unique name from metadata otherwise file path is used.
#' @param fileExtension \code{character}, file extension of mzXML formatted
#'  files. The directory is only searched for \emph{fileExtension} files.
#'  In most cases it would be \dQuote{"mzXML"} but sometimes you have to use
#'  \dQuote{xml}.
#' @param verbose \code{logical}, verbose output?
#'
#' @return A list of spectra.
#' \itemize{
#'  \item{\code{[[1]]spectrum$mass}: }{A vector of calculated mass.}
#'  \item{\code{[[1]]spectrum$intensity}: }{A vector of intensity values.}
#'  \item{\code{[[1]]metaData}: }{A list of metaData depending on read spectrum.}
#' }
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso \code{\link{readMzXmlFile}},
#' \code{\link[MALDIquantForeign]{importMzXml}}
#' @keywords IO
#' @rdname readMzXmlDir
#' @export
#' @examples
#'
#' ## load library
#' library("readMzXmlData")
#'
#' ## get examples directory
#' exampleDirectory <- system.file("Examples", package="readMzXmlData")
#'
#' ## read example spectra
#' spec <- readMzXmlDir(exampleDirectory)
#'
#' ## plot spectra
#' plot(spec[[1]]$spectrum$mass, spec[[1]]$spectrum$intensity, type="n")
#'
#' l <- length(spec)
#' legendStr <- character(l)
#' for (i in seq(along=spec)) {
#'   lines(spec[[i]]$spectrum$mass, spec[[i]]$spectrum$intensity, type="l",
#'         col=rainbow(l)[i])
#'   legendStr[i] <- basename(spec[[i]]$metaData$file)
#' }
#'
#' ## draw legend
#' legend(x="topright", legend=legendStr, col=rainbow(l), lwd=1)
#'
readMzXmlDir <- function(mzXmlDir, removeCalibrationScans=TRUE,
  removeMetaData=FALSE, rewriteNames=TRUE, fileExtension="mzXML",
  verbose=FALSE) {
  if (verbose) {
    message("Look for spectra in ", sQuote(mzXmlDir), " ...")
  }

  if ((!file.exists(mzXmlDir)) || (!file.info(mzXmlDir)$isdir)) {
    stop("Directory ", sQuote(mzXmlDir), " doesn't exists or is no ",
         "directory!")
  }

  ## look for mzXML files (alphabetical sort)
  files <- list.files(path=mzXmlDir,
                      pattern=paste0("^.*\\.", fileExtension, "$"),
                      recursive=TRUE)

  ## remove calibrations scans?
  if (removeCalibrationScans) {
    calibrationScans <- grep(pattern="[Cc]alibration", x=files, value=TRUE)
    if (length(calibrationScans) > 0) {
      files <- setdiff(files, calibrationScans)
    }
  }

  if (length(files) <= 0) {
    stop("Directory doesn't contain any ", fileExtension, " file.")
  }

  ## generate "path/files"
  files <- sapply(files, function(x) {
    x <- file.path(mzXmlDir, x)
    return(x)
  })

  ## read mzXML files
  mzXmlData <- list()
  for (i in seq(along=files)) {
    mzXmlFile <- .readMzXmlFile(mzXmlFile=files[i],
                                removeMetaData=removeMetaData, verbose=verbose)
    for (j in seq(along=mzXmlFile)) {
      spectra <- list()
      spectra$spectra <- mzXmlFile[[j]]
      mzXmlData <- c(mzXmlData, spectra)
    }
  }

  if (!removeMetaData & rewriteNames) {
    ## rewrite names
    if (verbose) {
      message("rewrite names ...")
    }

    names(mzXmlData) <- paste0("s", 1:length(mzXmlData))
  }

  return(mzXmlData)
}

