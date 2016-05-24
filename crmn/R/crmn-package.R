#' Normalize metabolomics data using CCMN and other methods
#'
#'
#' \tabular{ll}{
#' Package: \tab crmn \cr
#' Type: \tab Package \cr
#' Developed since: \tab 2009-05-14 \cr
#' Depends: \tab Biobase, pcaMethods (>= 1.20.2), pls, methods \cr
#' License: \tab GPL (>=3) \cr
#' LazyLoad: \tab yes \cr
#' }
#'
#' A package implementing the 'Cross-contribution compensating
#' multiple standard normalization'. Can be used to
#' normalize metabolomics data. Do \code{openVignette("crmn")} to see
#' the manual.
#'
#' @name crmn
#' @aliases crmn
#' @docType package
#' @title CRMN
#' @import Biobase
#' @import methods
#' @import pcaMethods
#' @author Henning Redestig
NULL

#' Mixture dilution series
#'
#' Multi-component dilution series. GC-TOF/MS measurements by Miyako Kusano.
#' Input concentrations are known and given in the original publication.
#' 
#' @name mix
#' @aliases mix
#' @usage data(mix)
#' @docType data
#' @title Dilution mixture dataset.
#' @examples
#'  data(mix)
#'  fData(mix)
#'  exprs(mix)
#'  pData(mix)
#' @author Henning Redestig
NULL

