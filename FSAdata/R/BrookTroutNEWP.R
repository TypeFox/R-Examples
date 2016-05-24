#' @title Catches in removal events for Brook Trout in the Nashwaak Experimental Watersheds Project.
#' 
#' @description Catches in removal events for Brook Trout (\emph{Salvelinus fontinalis}) in two streams in the the Nashwaak Experimental Watersheds Project on multiple dates.
#' 
#' @name BrookTroutNEWP
#' 
#' @docType data
#' 
#' @format A data frame of 16 observations on the following 7 variables:
#' \describe{
#'   \item{stream}{Stream (\code{UNM}=Upper Narrows Mountain Brook and \code{Hay}=Hyaden Brook).} 
#'   \item{section}{Section of stream.  See source.}
#'   \item{date}{Data of collections.}
#'   \item{first}{Catch on the first removal pass.} 
#'   \item{second}{Catch on the second removal pass.} 
#'   \item{third}{Catch on the third removal pass.} 
#'   \item{fourth}{Catch on the fourth removal pass.} 
#' }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Population size
#'    \item Abundance
#'    \item Removal
#'  }
#' 
#' @concept Abundance 'Population Size' Removal
#' 
#' @source From Table 1 in Schnute, J.  1983.  A new approach to estimating populations by the removal method.  Canadian Journal of Fisheries and Aquatic Sciences, 40:2153-2169.
#' 
#' @seealso See \code{\link{BrookTroutNEWP1}} for these data AND the results from Schnute (1983).
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BrookTroutNEWP)
#' 
#' ## extract data for one stream, section, and date (e.g., 3rd row)
#' BrookTroutNEWP[3,]
#' 
NULL
