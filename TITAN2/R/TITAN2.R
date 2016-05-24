#' Threshold Indicator Taxa Analysis
#'
#' Uses indicator species scores across binary partitions of a
#' sample set to detect congruence in taxon-specific changes of
#' abundance and occurrence frequency along an environmental
#' gradient as evidence of an ecological community threshold.
#'
#' Relevant references include: Baker, ME and RS King. 2010. A new
#' method for detecting and interpreting biodiversity and ecological
#' community thresholds.
#'
#' Methods in Ecology and Evolution 1(1): 25:37. King, RS and ME
#' Baker. 2010. Considerations for identifying and interpreting
#' ecological community thresholds.#' Journal of the North American
#' Benthological Association 29(3):998-1008. Baker ME and RS King.
#' 2013. Of TITAN and straw men: an appeal for greater understanding
#' of community data. Freshwater Science 32(2):489-506.
#'
#' @import parallel
#' @docType package
#' @name TITAN2
#' @importFrom graphics axis box legend mtext par plot points
#'   polygon segments symbols
#' @importFrom stats approxfun median quantile runif sd
#' @importFrom utils read.table write.table
#' @aliases TITAN2 package-TITAN2
NULL
