#' popEpi
#'
#' @name popEpi
#' @docType package
#' @title popEpi: Functions for large-scale epidemiological analysis
#' @author Janne Pitkaniemi, Joonas Miettinen, Karri Seppa, Matti Rantanen
#' @description 
#' \pkg{popEpi} is built for the needs of registry-based (large-scale)
#' epidemiological analysis. This is in most part enabled by the 
#' efficient \pkg{data.table} package for handling and aggregating large data sets. 
#' 
#' \pkg{popEpi} currently supplies some utility functions such as \code{\link{splitMulti}}
#' and \code{\link{get.yrs}} for preparing large data sets for epidemiological analysis.
#' Included are also a a few functions that can be used in 
#' epidemiological analysis such as \code{\link{sir}} for estimating
#' standardized incidence/mortality ratios (SIRs/SMRs) and \code{\link{survtab}} for 
#' estimating observed and relative/net survival as well as cumulative incidence
#' functions (CIFs).
#' 
#' Since there are many benefits to using \code{data.tables}, \pkg{popEpi} returns
#' outputs by default in the \code{data.table} format where appropriate. 
#' Since \code{data.table}
#' objects are usually modified by reference, this may have surprising side 
#' effects for users uninitiated in using \code{data.table}. To ensure
#' that appropriate outputs are in the \code{data.frame} format, set
#' \code{options("popEpi.datatable" = FALSE)}. However, \code{data.table}
#' usage is recommended due to better performance and testing coverage. 
#' \code{data.table} is used
#' by most functions internally in both cases.
#' 
#' 
NULL
