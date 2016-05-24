### Virtual classes
setClass("rawprof", representation="VIRTUAL")
setClass("fpmpi", contains="rawprof", representation="VIRTUAL")
setClass("mpip", contains="rawprof", representation="VIRTUAL")
setClass("tau", contains="rawprof", representation="VIRTUAL")


### Profiler class checker
valid_prof <- function(object)
{
  profilers <- c("fpmpi", "mpip", "tau")
  
  if ( !(object@profiler %in% profilers) )
    return("Invalid profiler")
  
  if (object@profiler != class(object@raw))
    return("'profiler' slot does not match class of 'raw' slot")
  
  return( TRUE )
}


#' Class prof
#' 
#' Class for Profiler Output
#' 
#' @slot profiler
#' The type of profiler used (e.g. fpmpi, mpiP). Stored as a character
#' vector.
#' @slot raw
#' The raw (non-processed) profiler output.  Storage is basically a
#' character vector, but set as a virtual class \code{rawprof}.
#' @slot parsed
#' A dataframe containing the processed version of the raw data.
#' 
#' @seealso \code{\link{read.prof}}
#' 
#' @aliases prof prof-class fpmpi-class mpip-class tau-class rawprof-class
#' 
#' @import methods
#' 
#' @name prof-class
#' @keywords Classes
#' @docType class
#' @exportClass prof
setClass(
  "prof", 
  representation(
    profiler="character",
    raw="rawprof",
    parsed="list"
  ),
  validity=valid_prof
)

