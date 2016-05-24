#' @title Ages of Morwong assigned from otoliths by two readers.
#' 
#' @description Ages assigned by two different readers to the otoliths of Jackass Morwong (\emph{Nemadactylus macropterus}).
#' 
#' @name Morwong3
#' 
#' @docType data
#' 
#' @format A data frame with 58 paired observations on the following 2 variables.
#'  \describe{
#'    \item{readerA}{Ages assigned by Reader A}
#'    \item{readerB}{Ages assigned by Reader B}
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
#' @source From Table 7  of Restall, J.E., and K. Krusic-Golub.  2004.  Development of jackass morwong age-length keys for 2000-2002. Final report to Australian Fisheries Management Authority. 13 pp. Primary Industries Research Victoria, Queenscliff.  Available at http://web-test.afma.gov.au/wp-content/uploads/2010/07/r03_1724b.pdf.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Morwong3)
#' str(Morwong3)
#' head(Morwong3)
#' plot(readerB~readerA,data=Morwong3)
#' with(Morwong3,table(readerA,readerB))
#' 
NULL