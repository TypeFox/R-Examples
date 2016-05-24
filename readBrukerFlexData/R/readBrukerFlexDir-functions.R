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

#' Reads recursively mass spectrometry data in Bruker Daltonics XMASS format.
#'
#' This function leads recursively all mass spectrometry data in
#' Bruker Daltonics XMASS format in a specified directory.
#'
#' @details
#' See \code{\link[readBrukerFlexData]{readBrukerFlexFile}}.
#'
#' @param brukerFlexDir \code{character}, path to \emph{directory} which
#'  should be read recursively.
#' @param removeCalibrationScans \code{logical}, if \code{TRUE} all scans in
#'  directories called \code{[Cc]alibration} will be ignored.
#' @param removeMetaData \code{logical}, to calculate mass data a lot of
#'  meta data are needed. To save memory they could be deleted after
#'  calculation.
#' @param useHpc \code{logical}, should Bruker Daltonics' High Precision
#'  Calibration be used if available? (see also:
#'  \code{\link[readBrukerFlexData]{.hpc}})
#' @param useSpectraNames \code{logical}, if \code{TRUE} all list elements
#'  get an unique name from metaData otherwise file path is used.
#'  (If \sQuote{removeMetaData} is \code{TRUE} \sQuote{useSpectraNames}
#'  has no effect.)
#' @param filterZeroIntensities \code{logical}, don't change it. If \code{TRUE}
#'  all intensities equal \code{0.0} are removed.
#'  (see also: \code{\link[readBrukerFlexData]{readBrukerFlexFile}})
#' @param verbose \code{logical}, print verbose messages?
#'
#' @return
#' A \code{list} of spectra.
#' \itemize{
#'     \item{\code{[[1]]$spectrum$mass}: }{A vector of calculated mass.}
#'     \item{\code{[[1]]$spectrum$intensity}: }{A vector of intensity values.}
#'     \item{\code{[[1]]$metaData}: }{A list of metaData depending on read spectrum.}
#' }
#'
#' @export
#' @seealso
#'  \code{\link[MALDIquantForeign]{importBrukerFlex}},
#'  \code{\link[readBrukerFlexData]{readBrukerFlexFile}},
#'  \code{\link[readBrukerFlexData]{.hpc}}
#' @rdname readBrukerFlexDir
#' @keywords IO
#' @examples
#' ## load library
#' library("readBrukerFlexData")
#'
#' ## get examples directory
#' exampleDirectory <- system.file("Examples", package="readBrukerFlexData")
#'
#' ## read example spectra
#' spec <- readBrukerFlexDir(file.path(exampleDirectory,
#'   "2010_05_19_Gibb_C8_A1"))
#'
#' ## plot spectra
#' plot(spec[[1]]$spectrum$mass, spec[[1]]$spectrum$intensity, type="n")
#'
#' l <- length(spec)
#' legendStr <- character(l)
#' for (i in seq(along=spec)) {
#'   lines(spec[[i]]$spectrum$mass, spec[[i]]$spectrum$intensity, type="l",
#'         col=rainbow(l)[i])
#'   legendStr[i] <- spec[[i]]$metaData$fullName
#' }
#'
#' ## draw legend
#' legend(x="topright", legend=legendStr, col=rainbow(l), lwd=1)
#'
readBrukerFlexDir <- function(brukerFlexDir, removeCalibrationScans=TRUE,
                              removeMetaData=FALSE, useHpc=TRUE,
                              useSpectraNames=TRUE,
                              filterZeroIntensities=FALSE, verbose=FALSE) {
  if (verbose) {
    message("Look for spectra in ", sQuote(brukerFlexDir), " ...")
  }

  if ((!file.exists(brukerFlexDir)) || (!file.info(brukerFlexDir)$isdir)) {
    stop("Directory ", sQuote(brukerFlexDir),
         " doesn't exists or is no directory!")
  }

  ## look for fid files (alphabetical sort)
  files <- list.files(path=brukerFlexDir, pattern="^fid$", recursive=TRUE)

  ## remove calibrations scans?
  if (removeCalibrationScans) {
    calibrationScans <- grep(pattern="[Cc]alibration", x=files, value=TRUE)
    if (length(calibrationScans) > 0) {
      files <- setdiff(files, calibrationScans)
    }
  }

  if (length(files) <= 0) {
    stop("Directory doesn't contain any fid file.")
  }

  ## generate 'path/files'
  files <- sapply(files, function(x) {
    x <- file.path(brukerFlexDir, x)
    return(x)
  })

  ## read fid files
  brukerFlexData <- lapply(X=files, FUN=function(f) {
    return(readBrukerFlexFile(fidFile=f, removeMetaData=removeMetaData,
                              useHpc=useHpc,
                              filterZeroIntensities=filterZeroIntensities,
                              verbose=verbose))
  })

  if (!removeMetaData & useSpectraNames) {
    ## rewrite names if metadata exists
    if (verbose) {
      message("look for spectra names ...")
    }

    names(brukerFlexData) <- sapply(X=brukerFlexData, FUN=function(x) {
      if (!is.null(x$metaData$sampleName)) {
        if (!is.na(x$metaData$sampleName)) {
          return(paste0("s", x$metaData$fullName, ".",
                        x$metaData$targetIdString))
        }
      }
      return(NA)
    }, USE.NAMES=FALSE)
  }

  return(brukerFlexData)
}

