#' Large deduplication data set
#'
#' Realization of RLBigData for deduplication of a single data set. Records are
#' stored as rows in \code{data}. Two records \code{data[i,]} and {data[j,]} are 
#' considered equal if and only if \code{identity[i]==identity[j]}
#'
#' @slot data Records to deduplicate
#' @slot identity Identity vector. 
setClass(
  Class = "RLResult",
  representation = representation(
    data = "RLBigData",
    prediction = "ff_vector"
  )#,
#  prototype = prototype(
#    data = NULL,
#    links = matrix(numeric(0), ncol=2, nrow=0),
#    possibleLinks = matrix(numeric(0), ncol=2, nrow=0),
#    nPairs = numeric(0)
#  )
)

# no constructor, is created by classifying methods

#setMethod(
#  f = "show",
#  signature = "RLResult",
#  definition = function(object)
#  {
#        
#  }
#)
#

