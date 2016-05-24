#                             edfReader package
#
# Purpose   :   Reading .edf(+)/.bdf(+) files
#
# Copyright :   (C) 2015-2016, Vis Consultancy, the Netherlands
#               This program is free software: you can redistribute it and/or modify
#               it under the terms of the GNU General Public License as published by
#               the Free Software Foundation, either version 3 of the License, or
#               (at your option) any later version.
#
#               This program is distributed in the hope that it will be useful,
#               but WITHOUT ANY WARRANTY; without even the implied warranty of
#               MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#               GNU General Public License for more details.
#
#               You should have received a copy of the GNU General Public License
#               along with edfReader package for R.  If not, see <http://www.gnu.org/licenses/>.
#
# History    :
#   Feb16 - Created, version 1.0.0
#   Mar16 - Version 1.1.0
#   Apr16 - Version 1.1.1, no changes
#
#' edfReader: A package for reading EDF(+) and BDF(+) files
#'
#' The edfReader package reads EDF(+) and BDF(+) files in two steps: first the header is read
#' and then the signals (using the header object as an parameter).
#'
#' @section edfReader functions:
#' \tabular{lll}{
#'   \code{\link{readEdfHeader}}  \tab \verb{ } \tab to read the file header with basic info about the signals \cr
#'   \code{\link{readEdfSignals}} \tab    \tab to read one or more recorded signals
#' }
#' The objects returned by these functions are described in the package vignette.
#'
#' @section Details:
#'  \tabular{lll}{
#'   Package \tab \verb{ } \tab edfReader \cr
#'   Version \tab  \tab 1.1.1 \cr
#'   Date \tab  \tab April 17, 2016 \cr
#'   Licence \tab \tab GPL version 3 or newer \cr
#'   GitHub \tab  \tab https://github.com/Pisca46/edfReader \cr
#'   Author \tab  \tab Jan Vis, Vis Consultancy \cr
#'   E-mail \tab  \tab jan@visconsultancy.eu \cr
#'   Web \tab  \tab visconsultancy.eu \cr
#' }
#' @section Acknowledgement:
#'    This package has used code from:
#'    \itemize{
#'      \item edf.R version 0.3 (27-11-2013) from Fabien Feschet, http://data-auvergne.fr/cloud/index.php/s/WYmFEDZylFWJzNs
#'      \item the work of Henelius Andreas as of July 2015, https://github.com/bwrc/edf
#'    }
#' @seealso
#'    For the vignette use the console command:\cr
#'    \code{vignette('edfReaderVignette', package = "edfReader")}\cr
#'    or click on \code{Index} below.
#
#' @aliases bdfReader
#' @docType package
#' @name edfReader
NULL

