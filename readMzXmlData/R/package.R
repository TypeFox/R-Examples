## Copyright 2012-2014 Sebastian Gibb
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

#' The readMzXmlData Package
#'
#' The package reads mass spectrometry data in mzXML format. \cr
#'
#' \tabular{ll}{
#'
#'  Package: \tab readMzXmlData\cr
#'
#'  Type: \tab Package\cr
#'
#'  Version: \tab 2.8.1\cr
#'
#'  Date: \tab 2015-09-16\cr
#'
#'  License: \tab GPL(version 3 or later)\cr
#'
#' }
#'
#' Main functions:
#'
#' \code{\link{readMzXmlFile}}: Reads mass spectrometry data in mzXML format.
#'
#' \code{\link{readMzXmlDir}}: Reads recursively mass spectrometry data in mzXML
#' format in a specific directory.
#'
#' \code{\link{mqReadMzXml}}: Reads mass spectrometry data into MALDIquant.
#'
#' @name readMzXmlData-package
#' @docType package
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso \code{\link{readMzXmlDir}}, \code{\link{readMzXmlFile}}
#' @references See website: \url{http://strimmerlab.org/software/maldiquant/}
#' @keywords IO
#'
#' @importFrom base64enc base64decode
#' @importFrom digest digest
#' @importFrom utils tail
#' @importFrom XML xmlEventParse
NULL

#' These functions are defunct and no longer available.
#'
#' \describe{
#'  \item{mqReadMzXml:}{use
#'    \code{\link[MALDIquantForeign]{importMzXml}} instead.}
#' }
#'
#' @title Removed functions in package \pkg{readMzXmlData}
#' @keywords internal
#' @name readMzXmlData-defunct
#' @aliases mqReadMzXml
#' @rdname readMzXmlData-defunct
#'
NULL

