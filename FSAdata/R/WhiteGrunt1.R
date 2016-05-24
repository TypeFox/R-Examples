#' @title Catch-at-age for White Grunt.
#' 
#' @description Catch-at-age for White Grunt (\emph{Haemulon plumierii}) collected from the central coast of Brazil.
#' 
#' @name WhiteGrunt1
#'
#' @docType data
#' 
#' @format A data frame with 25 observations on the following 2 variables.
#'  \describe{
#'    \item{age}{A numeric vector of assigned ages (from otoliths).}
#'    \item{catch}{A numeric vector of number of fish.} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Mortality 
#'    \item Catch curve 
#'  }
#'  
#' @concept Mortality 'Catch Curve'
#' 
#' @source From Figure 7 of Araujo, J.N. and A.S. Martins.  2007.  Age, growth and mortality of white grunt (\emph{Haemulon plumierii}) from the central coast of Brazil.  Scientia Marina 71:793-800.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WhiteGrunt1)
#' str(WhiteGrunt1)
#' head(WhiteGrunt1)
#' plot(log(catch)~age,data=WhiteGrunt1)
#' 
NULL