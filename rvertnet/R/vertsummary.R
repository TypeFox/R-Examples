#' Summarize a set of records downloaded from VertNet.
#'
#' Creates a simple summary of data returned by a VertNet search.
#'
#' @export
#' @param input Output from \code{\link{vertsearch}}, 
#'    \code{\link{searchbyterm}}, or \code{\link{spatialsearch}}. Required.
#' @param verbose Print progress and information messages. Default: TRUE
#' @return A list of summary statistics
#' @details \code{\link{vertsummary}} provides information on the sources, types and extent
#'    of data returned by a VertNet search.
#' @examples \dontrun{
#' recs <- vertsearch("Junco hyemalis")  # get occurrence records
#' vertsummary(recs)            # summarize occurrence records
#' 
#' vertsummary(vertsearch("Oncorhynchus clarki henshawi"))
#' }

vertsummary <- function(input, verbose = TRUE) vertsumwrapper(input, verbose)
