#' @title Catches in removal events for Brook Trout in the Nashwaak Experimental Watersheds Project.
#' 
#' @description Catches in removal events for Brook Trout (\emph{Salvelinus fontinalis}) in two streams in the the Nashwaak Experimental Watersheds Project on multiple dates.  Includes results from Schnute (1983).
#' 
#' @name BrookTroutNEWP1
#' 
#' @docType data
#' 
#' @format A data frame of 16 observations on the following 7 variables:
#' \describe{
#'   \item{sample}{A unique identified for the sample.}
#'   \item{stream}{Stream (\code{UNM}=Upper Narrows Mountain Brook and \code{Hay}=Hyaden Brook).} 
#'   \item{section}{Section of stream.  See source.}
#'   \item{date}{Data of collections.}
#'   \item{first}{Catch on the first removal pass.} 
#'   \item{second}{Catch on the second removal pass.} 
#'   \item{third}{Catch on the third removal pass.} 
#'   \item{fourth}{Catch on the fourth removal pass.}
#'   \item{Moran.N}{Schnute (1983) estimate of N using the Moran (1951) method.}
#'   \item{Moran.NLCI}{Schnute (1983) estimate of N 95\% LCI using the Moran (1951) method.}
#'   \item{Moran.NUCI}{Schnute (1983) estimate of N 95\% UCI using the Moran (1951) method.}
#'   \item{Moran.p}{Schnute (1983) estimate of p using the Moran (1951) method.}
#'   \item{Moran.LH}{Schnute (1983) negative log likelihood using the Moran (1951) method.}
#'   \item{Schnute.N}{Schnute (1983) estimate of N.}
#'   \item{Schnute.NLCI}{Schnute (1983) estimate of N 95\% LCI.}
#'   \item{Schnute.NUCI}{Schnute (1983) estimate of N 95\% UCI.}
#'   \item{Schnute.p1}{Schnute (1983) estimate of p1.}
#'   \item{Schnute.p}{Schnute (1983) estimate of p.}
#'   \item{Schnute.LH}{Schnute (1983) negative log-likelihood.}
#'   \item{ChiSq}{Schnute (1983) chi-square from likelihood ratio comparison of Moran and Schnute methods.}
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
#' @source From Tables 1-3 in Schnute, J.  1983.  A new approach to estimating populations by the removal method.  Canadian Journal of Fisheries and Aquatic Sciences, 40:2153-2169.
#' 
#' @seealso See \code{\link{BrookTroutNEWP}} for only the data (note the results from Schnute (1983)).
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BrookTroutNEWP1)
#' 
#' ## extract data for one stream, section, and date (e.g., 3rd row)
#' BrookTroutNEWP1[3,]
#' 
NULL
