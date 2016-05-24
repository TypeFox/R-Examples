#' @title Ages of American Shad assigned from scales by three readers at two times.
#'
#' @description Ages of American Shad (\emph{Alosa sapidissima}) assigned from scales by three readers at two times.
#' 
#' @details The true ages for fish in their sample were known because the Shad had been marked prior to being stocked. Additionally, 13 biologists twice (independently) estimated the age from scales for each fish.  The known age of the fish (\code{trueAge}) and the age estimates from three of the 13 biologists are available in this data.frame.  The estimated age variables are labeled with \code{ager}, a letter for the three biologists (\code{A}, \code{B}, or \code{C}) and a number for which time the scale was interpreted (\code{1} or \code{2}).  Some biologists chose not to assign an age to some scales and, thus, those data are missing (shown as \code{NA} values).
#'
#' @name ShadCR
#'
#' @docType data
#' 
#' @format A data frame with 53 observations on the following 8 variables.
#' \describe{
#'   \item{fishID}{A unique fish identification number}
#'   \item{trueAge}{The true age of the fish}
#'   \item{agerA1}{Ages assigned by reader A at time 1}
#'   \item{agerA2}{Ages assigned by reader A at time 2}
#'   \item{agerB1}{Ages assigned by reader B at time 1}
#'   \item{agerB2}{Ages assigned by reader B at time 2}
#'   \item{agerC1}{Ages assigned by reader C at time 1}
#'   \item{agerC2}{Ages assigned by reader C at time 2}
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
#' @note Used in the \href{http://derekogle.com/IFAR/}{Introductory Fisheries Analyses with R} book.
#' 
#' @source From McBride, R.S., Hendricks, M.L., and Olney, J.E. 2005. Testing the validity of Cating's (1953) method for age determination of American Shad using scales. Fisheries, 30:10-18.  Obtained directly from Rich McBride.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(ShadCR)
#' str(ShadCR)
#' head(ShadCR)
#' op <- par(mfrow=c(2,2),pch=19)
#' plot(agerA1~agerA2,data=ShadCR)
#' plot(agerB1~agerB2,data=ShadCR)
#' plot(agerC1~agerC2,data=ShadCR)
#' plot(agerA1~agerB1,data=ShadCR)
#' par(op)
NULL