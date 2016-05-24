#' Ages of Striped Bass assigned from scales and otoliths.
#' 
#' Ages of Striped Bass (\emph{Morone saxatilis}) assessed from heat-pressed scales and cracked-and-burnt otoliths.
#' 
#' @name StripedBass1
#' 
#' @docType data
#' 
#' @format A data frame of 343 observations on the following 2 variables:
#'  \describe{
#'    \item{ageO}{Age assigned from examinaton of otoliths}
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
#' @concept Age Precision Bias
#' 
#' @source From Figure 6 in Chapter 10 (Striped Bass) of the VMRC Final Report on Finfish Ageing, 2000 by the Center for Quantitative Fisheries Ecology at Old Dominion University.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(StripedBass1)
#' str(StripedBass1)
#' head(StripedBass1)
#' plot(ageS~ageO,data=StripedBass1)
#' xtabs(~ageO+ageS,data=StripedBass1)
#' 
NULL
