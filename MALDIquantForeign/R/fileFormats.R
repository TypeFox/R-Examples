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

#' Supported file formats
#'
#' This function prints all file formats supported by
#' \code{\link{MALDIquantForeign-package}}.
#'
#' \subsection{Import}{
#'
#' \tabular{ll}{
#'  txt \tab \code{\link[MALDIquantForeign]{importTxt}} \cr
#'  tab \tab \code{\link[MALDIquantForeign]{importTab}} \cr
#'  csv \tab \code{\link[MALDIquantForeign]{importCsv}} \cr
#'  fid \tab \code{\link[MALDIquantForeign]{importBrukerFlex}} \cr
#'  ciphergen \tab \code{\link[MALDIquantForeign]{importCiphergenXml}} \cr
#'  mzXML \tab \code{\link[MALDIquantForeign]{importMzXml}} \cr
#'  mzML \tab \code{\link[MALDIquantForeign]{importMzMl}} \cr
#'  imzML \tab \code{\link[MALDIquantForeign]{importImzMl}} \cr
#'  analyze \tab \code{\link[MALDIquantForeign]{importAnalyze}} \cr
#'  cdf \tab \code{\link[MALDIquantForeign]{importCdf}} \cr
#'  msd \tab \code{\link[MALDIquantForeign]{importMsd}} \cr
#' }
#' }
#'
#' \subsection{Export}{
#'
#' \tabular{ll}{
#'  tab \tab \code{\link[MALDIquantForeign]{exportTab}} \cr
#'  csv \tab \code{\link[MALDIquantForeign]{exportCsv}} \cr
#'  imzML \tab \code{\link[MALDIquantForeign]{exportImzMl}} \cr
#'  msd \tab \code{\link[MALDIquantForeign]{exportMsd}} \cr
#'  mzML \tab \code{\link[MALDIquantForeign]{exportMzMl}} \cr
#' }
#' }
#'
#' @return a \code{list} with two named elements (\code{import} and
#' \code{export}) containing a \code{character} vector of supported file types.
#'
#' @seealso
#'  \code{\link[MALDIquantForeign]{export}},
#'  \code{\link[MALDIquantForeign]{import}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}
#' @examples
#' library("MALDIquantForeign")
#'
#' supportedFileFormats()
#'
#' @rdname supportedFileFormats-functions
#' @export
supportedFileFormats <- function() {
  list(import=importFormats$type,
       export=exportFormats$type)
}

importFormats <- data.frame(type=c("txt", "tab", "csv", "fid", "ciphergen",
                                   "mzxml", "mzml", "imzml", "analyze", "cdf",
                                   "msd"),
                            pattern=c("^.*\\.txt$", "^.*\\.tab$",
                                      "^.*\\.csv$", "^fid$",
                                      "^.*\\.xml$", "^.*\\.mzXML$",
                                      "^.*\\.mzML$", "^.*\\.imzML$",
                                      "^.*\\.hdr$", "^.*\\.cdf$", "^.*\\.msd$"),
                            handler=c(rep(".importTab", 2),
                                      ".importCsv", ".importBrukerFlex",
                                      ".importCiphergenXml", ".importMzXml",
                                      ".importMzMl", ".importImzMl",
                                      ".importAnalyze", ".importCdf",
                                      ".importMsd"),
                            stringsAsFactors=FALSE)

exportFormats <- data.frame(type=c("tab", "csv", "msd", "mzml", "imzml"),
                            extension=c("tab", "csv", "msd", "mzML", "imzML"),
                            onefile=c(FALSE, FALSE, FALSE, TRUE, TRUE),
                            handler=c(".exportTab", ".exportCsv",
                                      ".exportMsd", ".exportMzMl",
                                      ".exportImzMl"),
                            stringsAsFactors=FALSE)
