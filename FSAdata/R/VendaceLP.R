#' @title Stock and recruitment data for Vendace from Lake Puulavesi, 1982-1996.
#' 
#' @description Vendace (\emph{Coregonus albula}) recruitment by year in Lake Puulavesi, 1982-1996.
#' 
#' @name VendaceLP
#' 
#' @docType data
#' 
#' @format A data frame of 15 observations on the following 3 variables:
#'  \describe{
#'    \item{year}{Year of data} 
#'    \item{stock}{Spawning stock index -- autumn biomass (kg/ha) of age-1+ fish}
#'    \item{recruits}{Recuit index -- density (indivs/ha) of age-0+ fish in first autumn}
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
#' @source From (approximately) Figure 1 and 2 of Marjomaki, T.J. 2004.  Analysis of the spawning stock-recruitment relationship of vendace (\emph{Coregonus albula} (L.)) with evaluation of alternative models, additional variables, biases and errors.  Ecology of Freshwater Fish 13:46-60.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(VendaceLP)
#' str(VendaceLP)
#' head(VendaceLP)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(recruits~year,data=VendaceLP,type="l")
#' plot(recruits~stock,data=VendaceLP)
#' par(op)
#' 
NULL