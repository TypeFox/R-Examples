#' @title Annual catches of yearling Walleye in bottom trawls from Lake Winnebago, WI, 1986-2010.
#' 
#' @description Annual catches of yearling Walleye (\emph{Sander vitreus}) in bottom trawls from Lake Winnebago, WI, 1986-2010.
#' 
#' @details The catch of yearling Walleye and number of trawl tows by year are in this data.frame. The CPE is catch divided by number of tows. Koenigs et al. (2015) rescaled the CPE values to have a mean of 0 and a standard deviation of 1.

#' @name WalleyeWyrlng
#' 
#' @docType data
#' 
#' @format A data frame with 35 observations on the following 4 variables.
#'  \describe{
#'    \item{tows}{Number of trawl tows (i.e., effort)}
#'    \item{year}{Year of capture}
#'    \item{yearlings}{Number of yearling Walleye captured}
#'    \item{yrclass}{Year-class of the captured yearlings (capture year minus 1)}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Year-class Strength
#'    \item Recruitment
#'  }
#'  
#' @concept 'Year-class Strength' Recruitment
#' 
#' @note Used in the \href{http://derekogle.com/IFAR/}{Introductory Fisheries Analyses with R} book.
#' 
#' @source From Koenigs, R.P., Bruch, R.M., Stelzer, R.S., and Kamke, K.K. 2015. Validation of otolith ages for Walleye (\emph{Sander vitreus}) in the Winnebago System.  Fisheries Research, 167:13-21.  Obtained directly from Ryan Koenigs.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WalleyeWyrlng)
#' str(WalleyeWyrlng)
#' head(WalleyeWyrlng)
#' plot(yearlings~yrclass,data=WalleyeWyrlng)
#' 
NULL