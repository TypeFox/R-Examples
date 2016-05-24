#' Class \code{"rtadpcr"} - real-time array digital PCR experiments
#' 
#' A class designed to contain results from real-time array digital PCR
#' experiments. Data is represented as matrix, where each column describes
#' different measurement point (i.e. cycle number) and every row different
#' partition.
#' 
#' 
#' @name rtadpcr-class
#' @aliases rtadpcr-class rtadpcr
#' @docType class
#' @section Slots: \describe{ \item{list(".Data")}{\code{"matrix"} containing
#' data from array. See Description.}\item{:}{\code{"matrix"} containing data
#' from array. See Description.} \item{list("n")}{Object of class
#' \code{"integer"} equal to the number of partitions.}\item{:}{Object of class
#' \code{"integer"} equal to the number of partitions.}
#' \item{list("type")}{Object of class \code{"character"} defining type of
#' data.}\item{:}{Object of class \code{"character"} defining type of data.} }
#' @author Michal Burdukiewicz.
#' @seealso End-point array digital PCR: \code{\linkS4class{adpcr}}.
#' 
#' Droplet digital PCR: \code{\linkS4class{ddpcr}}.
#' @keywords classes real-time
#' @examples
#' 
#' #none
#' 
NULL

setClass("rtadpcr", contains = "matrix", representation(.Data = "matrix",
                                                      n = "integer",
                                                      type = "character"))