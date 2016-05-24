#' @title Ages of Atlantic Croaker assigned from otoliths by two readers.
#' 
#' @description Otolith age of Atlantic croaker (\emph{Micropogonias undulatus}) from two readers.
#' 
#' @name Croaker1
#' 
#' @docType data
#' 
#' @format A data frame of 317 observations on the following 2 variables:
#'  \describe{
#'    \item{reader1}{Age assigned by the first reader}
#'    \item{reader2}{Age assigned by the second reader} 
#'  }
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
#' @source From Figure 2 in Chapter 1 (Atlantic Croaker) of the VMRC Final Report on Finfish Ageing, 1999 by the Center for Quantitative Fisheries Ecology at Old Dominion University.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Croaker1)
#' str(Croaker1)
#' head(Croaker1)
#' plot(reader2~reader1,data=Croaker1)
#' xtabs(~reader1+reader2,data=Croaker1)
#' 
NULL