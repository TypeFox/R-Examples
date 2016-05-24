#' @title Catches in removal events of Cutthroat Trout and Coho Salmon in Little Stawamus Creek (British Columbia, Canada) in 1997.
#' 
#' @description Catches of Coho Salmon (\emph{Oncorhynchus kisutch}) and Cutthroat Trout (\emph{Oncorhynchus clarki}) in consecutive removal events at various locations in Little Stawamus Creek (British Columbia, Canada) in 1997.
#' 
#' @name Deckeretal1999
#' 
#' @docType data
#' 
#' @format A data frame of 26 observations on the following 10 variables:
#'  \describe{
#'    \item{reach}{Reach number of sampling location.}
#'    \item{habitat}{Habitat type of sampling location -- \code{pool}, \code{riffle}, or \code{run}.}
#'    \item{strata}{Stratum number of sampling location.}
#'    \item{area}{Area (m^2) of sampling location.}
#'    \item{coho1}{Coho Salmon removed on the first pass.}
#'    \item{coho2}{Coho Salmon removed on the second pass.}
#'    \item{coho3}{Coho Salmon removed on the third pass.}
#'    \item{cutt1}{Cutthroat Trout removed on the first pass.}
#'    \item{cutt2}{Cutthroat Trout removed on the second pass.}
#'    \item{cutt3}{Cutthroat Trout removed on the third pass.}
#'  }
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
#' @source From Appendix 2a and 2b in Decker, A.S., J.M. Bratty, S.C. Riley, and J. Korman.  1999.  Estimating standing stock of juvenile coho salmon (\emph{Oncorhynchus kisutch}) and cutthroat trout (\emph{Oncorhynchus clarki}) in a small stream: a comparison of sampling designs.  Canadian Technical Report of Fisheries and Aquatic Sciences 2282.  24 pp.  [Was (is?) from http://www.dfo-mpo.gc.ca/Library/239481.pdf.]
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Deckeretal1999)
#' str(Deckeretal1999)
#' head(Deckeretal1999)
#' 
#' ## extract data for one sampling location (e.g., 3rd row)
#' Deckeretal1999[3,]
#' 
NULL