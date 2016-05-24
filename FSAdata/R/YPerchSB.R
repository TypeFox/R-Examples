#' @title Stock and recruitment data for Yellow Perch from South Bay, Lake Huron, 1950-1983.
#' 
#' @description Yellow Perch (\emph{Perca flavescens}) stock and recruitment by year in South Bay, Lake Huron, 1950-1983.
#' 
#' @name YPerchSB
#' 
#' @docType data
#' 
#' @format A data frame with 15 observations on the following 3 variables:
#'  \describe{
#'    \item{year}{Year of data}
#'    \item{stock}{Spawning stock (number per set)}
#'    \item{recruits}{Recruits (number per set)}
#'  }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Stock-Recruit
#'    \item Recruitment
#'  }
#' 
#' @concept 'Stock-Recruit' Recruitment
#' 
#' @source From Henderson, B.A. 1985. Factors affecting growth and recruitment of yellow perch, \emph{Perca flavescens} Mitchill, in South Bay, Lake Huron. Journal of Fisheries Biology 26:449-458.  Obtained from Ransom Myers online database which was (is?) at http://ram.biology.dal.ca/~myers/data.html.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(YPerchSB)
#' str(YPerchSB)
#' head(YPerchSB)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(recruits~year,data=YPerchSB,type="b")
#' plot(recruits~stock,data=YPerchSB)
#' par(op)
#' 
NULL