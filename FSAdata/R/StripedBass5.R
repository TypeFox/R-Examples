#' Ages of Striped Bass assigned from otoliths by two readers.
#' 
#' Ages assigned by two different readers to the otoliths of Striped Bass (\emph{Morone saxatilis}).
#' 
#' @name StripedBass5
#' 
#' @docType data
#' 
#' @format A data frame of 458 observations on the following 2 variables:
#'  \describe{
#'    \item{reader1}{Ages assigned by the first reader}
#'    \item{reader2}{Ages assigned by the second reader} 
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
#' @seealso \code{\link{StripedBass4}} and \code{\link{StripedBass6}}.
#' 
#' @source From Figure 5 in Chapter 10 (Striped Bass) of the VMRC Final Report on Finfish Ageing, 2003 by the Center for Quantitative Fisheries Ecology at Old Dominion University.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(StripedBass5)
#' str(StripedBass5)
#' head(StripedBass5)
#' plot(reader2~reader1,data=StripedBass5)
#' xtabs(~reader1+reader2,data=StripedBass5)
#' 
NULL