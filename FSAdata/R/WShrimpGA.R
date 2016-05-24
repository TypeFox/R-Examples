#' @title Stock and recruitment data for White Shrimp off the coast of Georgia (USA), 1979-2000.
#' 
#' @description White Shrimp (\emph{Litopenaeus setiferus}) stock and recruitment data from off the coast of Georgia (USA), 1979-2000.
#' 
#' @note No graph was shown for stock values vs. year so some stock values had to be haphazardly paired with recruit values -- especially for stock values between 2.5 and 3.0.
#' 
#' @name WShrimpGA
#' 
#' @docType data
#' 
#' @format A data frame with 22 observations on the following 3 variables:
#'  \describe{
#'    \item{year}{Year of data}
#'    \item{stock}{Spawning stock index -- CPUE in June assessment surveys}
#'    \item{recruits}{Recuit index -- commercial landings in pounds from Aug-Jan}
#'    \item{yrconf}{A code of confidence for whether or not the stock value is known to come from the year shown (see the note)}
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
#' @source From (approximately) Figures 2 and 3 of Belcher, C.N., and C.A. Jennings.  2004.  Evaluation of stock-recruitment curves for white shrimp in Georgia.  North American Journal of Fisheries Management. 24:654-661.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WShrimpGA)
#' str(WShrimpGA)
#' head(WShrimpGA)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(recruits~year,data=WShrimpGA,type="b")
#' plot(recruits~stock,data=WShrimpGA)
#' par(op)
#' 
NULL