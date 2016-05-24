#' Stock and recruitment data for Pacific Sardines, 1935-1990.
#'
#' Pacific Sardine (\emph{Sardinops sagax}) stock and recruitment by year, 1935-1990.
#'
#' @name SardinesPacific
#' @docType data
#' 
#' @format A data frame of 34 observations on the following 3 variables:
#' \describe{
#'   \item{year}{Year of data}
#'   \item{ssb}{Spawning stock biomass -- millions of fish}
#'   \item{recruits}{Recruitment index -- millions of fish}
#' }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Stock-Recruit
#'    \item Recruitment
#'  }
#' 
#' @concept 'Stock-Recruit' Recruitment
#' 
#' @source From Jacobson, L.D. and A.D. MacCall.  1995.  Stock-recruitment models for Pacific sardine (\emph{Sardinops sagax}). Canadian Journal of Fisheries and Aquatic Sciences. 52:566-577.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(SardinesPacific)
#' str(SardinesPacific)
#' head(SardinesPacific)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(recruits~year,data=SardinesPacific,type="l")
#' plot(recruits~ssb,data=SardinesPacific)
#' par(op)
#'
NULL