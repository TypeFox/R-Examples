#' @title Catches-at-age for male and female Walleye from Lake Winnebago, WI, 2010.
#' 
#' @description Catches-at-age for male and female Walleye from Lake Winnebago, WI, 2010.
#' 
#' @details Koenigs et al. (2015) captured adult Walleye from Lake Winnebago during spawning assessments in 2010. The sex was recorded and ages were estimated from sectioned otoliths for each fish. Koenigs et al. (2015) fit separate catch curves to female and male Walleye.

#' @name WalleyeWad
#' 
#' @docType data
#' 
#' @format A data frame with 18 observations on the following 3 variables.
#'  \describe{
#'    \item{age}{Age (yrs; from pectoral fin ray)}
#'    \item{numF}{Number of captured females}
#'    \item{numM}{Number of captured males}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Mortality
#'    \item Catch Curve
#'  }
#'  
#' @concept 'Catch Curve' Mortality
#' 
#' @note Used in the \href{http://derekogle.com/IFAR/}{Introductory Fisheries Analyses with R} book.
#' 
#' @source From Koenigs, R.P., Bruch, R.M., Stelzer, R.S., and Kamke, K.K. 2015. Validation of otolith ages for Walleye (\emph{Sander vitreus}) in the Winnebago System.  Fisheries Research, 167:13-21.  Obtained directly from Ryan Koenigs.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WalleyeWad)
#' str(WalleyeWad)
#' head(WalleyeWad)
#' plot(numF~age,data=WalleyeWad)
#' points(numM~age,data=WalleyeWad,pch=19)
#' 
NULL