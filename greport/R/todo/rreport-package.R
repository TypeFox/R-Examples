#' This package creates reports.
#'
#' @author Frank E Harrell Jr \email{f.harrell@@vanderbilt.edu}
#'
#' Maintainer: Charles Dupont \email{charles.dupont@@vanderbilt.edu}
#'
#' @importFrom chron as.chron dates seq.dates chron years days hours minutes seconds
#' @importFrom lattice bwplot
#' @importFrom rms Surv survfit.formula survplot.survfit
#' @import Hmisc
#' @docType package
#' @aliases rreport package-rreport
#' @name rreport
NULL

# The caching and check for conflicts require looking for a pattern of objects; the search may be avoided by defining an object ‘.noGenerics’
# see ?library
.noGenerics <- TRUE

#' RReport Package Options
#'
#' @aliases rreport.options
#' @section Options used in rreport:
#' \describe{
#'  \item{\code{rreport.gtype}:}{graphing device (ps, pdf, interactive)}
#'  \item{\code{rreport.appendix.file.name}:}{filename for appendix}
#'  \item{\code{rreport.generated.tex.dir}:}{directory name for report tex files}
#'  \item{\code{rreport.graphics.dir}:}{directory name for report graphic files}
#'  \item{\code{rreport.closed.filename.mask}:}{mask for closed report filenames}
#'  \item{\code{rreport.open.filename.mask}:}{mask for open report filenames}
#' }
#' @name pkgOptions
#' @seealso \code{\link[base]{options}}
NULL

.defaultRreportOptions <- function() {
  list(
    rreport.gtype = 'pdf',
    rreport.appendix.file.name = 'app.tex',
    rreport.generated.tex.dir = 'gentex',
    rreport.graphics.dir = 'pdf',
    rreport.closed.filename.mask = NULL,
    rreport.open.filename.mask = 'O%s'
  )
}

.onLoad <- function(libname, pkgname) {
  options(.defaultRreportOptions())
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("rreport library by Frank E Harrell Jr\n\nType library(help='rreport') to see overall documentation.\n\n")
}

# remove rreport options
.onUnload <- function(libpath) {
  ropts <- grep("^rreport", names(options()), value=TRUE)
  nulls <- vector('list', length(ropts))
  names(nulls) <- ropts
  options(nulls)
}

# questions?
# accrualReport: no visible binding for global variable ‘code.infig’
# aeReport: no visible binding for global variable ‘weeks’
# completenessReport: no visible binding for global variable ‘compFullCaptionDone’
# freqReport: no visible binding for global variable ‘name’
# rangeCheck: no visible binding for global variable ‘dataframe’
# there is no "gtype" (interactive/pdf/ps)
# endPlot
# putFig
# startPlot

# mixedvarReport: no visible global function definition for ‘Key’

# S3 methods shown with full name in documentation object 'getReferenceObject':
#   ‘print.latexReference’
# 
# S3 methods shown with full name in documentation object 'floor.chron':
#   ‘floor.chron’ ‘ceiling.chron’
# 
# The \usage entries for S3 methods should use the \method markup and not their full name.
# See the chapter ‘Writing R documentation files’ in the ‘Writing R Extensions’ manual.
# * checking Rd contents ... WARNING
# Argument items with no description in Rd object 'getReferenceObject':
#   ‘refD’ ‘newMarker’ ‘keyword’ ‘label’
#
# regarding survReport
#' This report assumes units are in days.
#' I can use this, but the output looks dumb.
#' valueUnit(mydata$time) <- "Month"
