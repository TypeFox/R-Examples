#' @title Ages of Muskellunge assigned from scales and cleithra.
#' 
#' @description Ages of St. Lawrence River, ONT, Muskellunge (\emph{Esox masquinongy}) assessed from scales and cleithra.
#' 
#' @name MuskieSLR
#' 
#' @docType data
#' 
#' @format A data frame of 43 observations on the following 2 variables:
#'  \describe{
#'    \item{ageC}{Age assigned from examinaton of cleithrum}
#'    \item{ageS}{Age assigned from examination of scales} 
#'  }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Age Comparison
#'    \item Age Precision 
#'    \item Age Bias
#'    \item Ageing Error
#'  }
#' 
#' @concept Age Precision Bias 'Age Comparison'
#' 
#' @source From Figure 3 in Casselman,J.M.  1983. Age and growth assessment of fish from their calcified structures - techniques and tools. In E.D. Prince and L.M. Pulos, editors, Proceedings of the international workshop on age determination of oceanic pelagic fishes: Tunas, billfishes, and sharks, volume NOAA Technical Report, NMFS 8:1-17.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(MuskieSLR)
#' str(MuskieSLR)
#' head(MuskieSLR)
#' plot(ageS~ageC,data=MuskieSLR)
#' xtabs(~ageC+ageS,data=MuskieSLR)
#' 
NULL