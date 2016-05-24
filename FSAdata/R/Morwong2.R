#' @title Ages of Morwong assigned from otoliths by Reader B at two times.
#' 
#' @description Ages assigned at two different times by Reader B to the otoliths of Jackass Morwong (\emph{Nemadactylus macropterus}).
#' 
#' @name Morwong2
#' 
#' @docType data
#' 
#' @format A data frame with 136 observations on the following 2 variables.
#'  \describe{
#'    \item{first}{Ages assigned on the first reading}
#'    \item{second}{Ages assigned on the second reading}
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
#' @source From Table 6 of Restall, J.E., and K. Krusic-Golub.  2004.  Development of jackass morwong age-length keys for 2000-2002. Final report to Australian Fisheries Management Authority. 13 pp. Primary Industries Research Victoria, Queenscliff.  Available at http://web-test.afma.gov.au/wp-content/uploads/2010/07/r03_1724b.pdf.  
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Morwong2)
#' str(Morwong2)
#' head(Morwong2)
#' plot(second~first,data=Morwong2)
#' xtabs(~first+second,data=Morwong2)
#' 
NULL