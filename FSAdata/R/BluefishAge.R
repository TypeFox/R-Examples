#' @title Ages of Bluefish assigned from otoliths by two readers.
#'
#' @description Ages assigned to Bluefish (\emph{Pomatomus saltatrix}) otoliths by two readers.
#'
#' @name BluefishAge
#'
#' @docType data
#' 
#' @format A data frame with 445 observations on the following 2 variables.
#' \describe{
#'   \item{r1}{Ages assigned by the first reader}
#'   \item{r2}{Ages assigned by the second reader}
#' }
#' @section Topic(s): \itemize{
#'   \item Age Comparison
#'   \item Age Precision 
#'   \item Age Bias
#'   \item Ageing Error
#'  }
#' 
#' @concept Age Precision Bias 'Age Comparison'
#' 
#' @source From Figure 2 in Chapter 3 (Bluefish) of the VMRC Final Report on Finfish Ageing, 2003 by the Center for Quantitative Fisheries Ecology at Old Dominion University.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BluefishAge)
#' str(BluefishAge)
#' head(BluefishAge)
#' plot(r1~r2,data=BluefishAge)
#'
NULL
