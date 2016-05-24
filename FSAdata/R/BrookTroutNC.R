#' @title Stock and recruitment data for Brook Trout from Ball Creek, NC, 1991-2004.
#' 
#' @description Stock and recruitment data for Brook Trout (\emph{Salvelinus fontinalis}) from Ball Creek, NC, 1991-2004.
#' 
#' @note The authors fit a linear model to the stock-recruit relationship.
#' 
#' @name BrookTroutNC
#' 
#' @docType data
#' 
#' @format A data frame with 10 observations on the following 2 variables.
#'  \describe{
#'    \item{adult}{a numeric vector giving autumn adult density (number per square meter)}
#'    \item{yoy}{a numeric vector giving autumn YOY density (number per square meter) in following year} 
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
#' @source From (approximately) Figure 5 in Grossman, G.D., R.E. Ratajczak, C.M. Wagner, and J.T. Petty.  2010. Dynamics and regulation of the southern brook trout (\emph{Salvelinus fontinalis}) population in an Appalachian stream.  Freshwater Biology 55:1494-1508.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BrookTroutNC)
#' str(BrookTroutNC)
#' head(BrookTroutNC)
#' plot(adult~yoy,data=BrookTroutNC)
#' 
NULL
