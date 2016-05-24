#' treeplyr: 'dplyr' Functionality for Matched Tree and Data Objects
#' 
#' Matches phylogenetic trees and trait data, and allows simultaneous 
#' manipulation of the tree and data using 'dplyr'.
#' 
#' @docType package
#' @author Josef Uyeda
#' @name treeplyr
#' @useDynLib treeplyr
#' @import ape dplyr Rcpp
#' @importFrom lazyeval all_dots
#' @importFrom phytools phylosig
#' @importFrom geiger fitContinuous
#' @importFrom grDevices rainbow
#' @importFrom graphics legend lines locator par plot points text
#' @importFrom stats setNames
#' @importFrom utils type.convert
NULL
