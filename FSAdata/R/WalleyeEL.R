#' @title Stock and recruitment data for Walleye from Escanaba Lake, WI, 1958-1992.
#' 
#' @description Abundance of age-0 and age-5 and older Walleye (\emph{Sander vitreus}), abundance of adult Yellow Perch (\emph{Perca flavescens}), and coefficient of variation of May temperatures for Escanaba Lake, WI, 1958-1992.
#' 
#' @name WalleyeEL
#' 
#' @docType data
#' 
#' @format A data frame of 39 observations on the following 5 variables:
#'  \describe{
#'    \item{yrclass}{Year-class of the data}
#'    \item{age0}{Abundance of age-0 Walleye (recruits)}
#'    \item{age5}{Abundance of age-5 and older Walleye (stock)}
#'    \item{maycv}{Coefficient of variation of May temperatures in birth year}
#'    \item{yep}{Abundance of adult (larger than 152.4 mm) Yellow Perch}
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
#' @source Hansen, M. J., M. A. Bozek, J. R. Newby, S. P. Newman, and M. D. Staggs.  1998. Factors affecting recruitment of walleyes in Escanaba Lake, Wisconsin, 1958-1995.  North American Journal of Fisheries Management 18:764-774.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WalleyeEL)
#' str(WalleyeEL)
#' head(WalleyeEL)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(age0~yrclass,data=WalleyeEL,type="l")
#' plot(age0~age5,data=WalleyeEL)
#' par(op)
#' 
NULL