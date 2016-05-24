#' @title Stock and recruitment data for Lake Huron Bloaters, 1981-1996.
#' 
#' @description Egg deposition and relative abundance of age-3 Lake Huron Bloaters (\emph{Coregonus hoyi}) by year, 1981-1996.
#' 
#' @name BloaterLH
#' 
#' @docType data
#' 
#' @format A data frame of 16 observations on the following 3 variables:
#'  \describe{
#'    \item{year}{Year of data (1981-1996)}
#'    \item{eggs}{Millions of eggs deposited} 
#'    \item{age3}{Relative abundance of age-3 fish}
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
#' @source From (approximately) Figure 7 of Schaeffer, J.S. 2004.  Population dynamics of bloaters \emph{Coregonus hoyi} in Lake Huron, 1980-1998.  Ann Zool Fennici. 41:271-279.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BloaterLH)
#' str(BloaterLH)
#' head(BloaterLH)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(eggs~year,data=BloaterLH,type="l")
#' plot(eggs~age3,data=BloaterLH)
#' par(op)
#' 
NULL
