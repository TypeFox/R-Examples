.onAttach <- function(libname, pkgname) {
   packageStartupMessage("This is hetmeta 0.1.0. For an overview type: help('hetmeta-package').")
}

#' Heterogeneity Measures In Meta-Analysis
#' @name hetmeta-package
#' @docType package
#' @import metafor
#' @importFrom stats qnorm
#'
#' @description
#' The \code{hetmeta} package contains functions useful to assess the presence and to quantifying
#' the impact of statistical heterogeneity.
#' Several measures of heterogeneity are implemented in the \code{\link{hetmeta}} function.
#'
#' All the functions in the packages requires a meta-analytic model of class \code{\link{rma.uni}}
#' that can be easily obtained using the \code{\link{metafor}} package.
#' See \code{\link{metafor-package}} for a comprehensive and detailed description.
#' @section Functions and data included in the package:
#'
#' @section Functions and data included in the package:
#'
#' The main function is \code{\link{hetmeta}}, which calculates the measures of heterogeneity
#' in an object of class "\code{hetmeta}" (see \code{\link{hetmetaObject}}). The methods
#' \code{\link{print.hetmeta}} and \code{\link{confint.hetmeta}} defines function for
#' printing results and deriving confidence intervals.
#'
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#'
#' @references
#'
#' Crippa A, Khudyakov P, Wang M, Orsini N, Spiegelman D. A new measure of between-studies heterogeneity in meta-analysis. 2016. \emph{Stat. Med.} In Press.
#'
#' DerSimonian R, Laird N. Meta-analysis in clinical trials. \emph{Control. Clin. Trials} 1986; 7(3):177-188.
#'
#' Rebecca HJ, Thompson J. Detecting and describing heterogeneity in meta-analysis.
#' \emph{Stat. Med.} 17.8 (1998): 841-856.
#'
#' Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. \emph{Stat. Med.} 2002; 21(11):1539-1558.

NULL

