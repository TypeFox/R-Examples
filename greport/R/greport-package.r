#' Graphical Reporting for Clinical Trials
#'
#' @author Frank E Harrell Jr \email{f.harrell@@vanderbilt.edu}
#'
#' @maintainer Frank E Harrell Jr \email{f.harrell@@vanderbilt.edu}
#'
#' @export Merge accrualReport dNeedle dReport eReport endPlot exReport getgreportOption nriskReport putFig sampleFrac setgreportOption startPlot survReport
#' @import Hmisc ggplot2 lattice data.table methods
#' @importFrom latticeExtra useOuterStrips
#' @importFrom rms npsurv survplot
#' @importFrom survival Surv survfit
#' @importFrom Formula Formula model.part
#' @importFrom grDevices adjustcolor dev.off gray pdf
#' @importFrom graphics abline axis box grconvertX grconvertY lines par plot plot.new points text
#' @importFrom stats as.formula median model.frame qnorm reshape sd terms
#' @docType package
#' @aliases greport package-greport
#' @name greport

# The caching and check for conflicts require looking for a pattern of objects; the search may be avoided by defining an object .noGenerics
# see ?library
.noGenerics <- TRUE
