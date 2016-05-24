#' @title Recruitment time-series for Rainbow Smelt in Lake Erie, 1977-1996.
#' 
#' @description Rainbow Smelt (\emph{Osmerus mordax}) recruitment time series for Lake Erie, 1977-1996.
#' 
#' @note Zeros were changed to 0.1 in 1984, 1987, 1991, and 1994.
#' 
#' @name RBSmeltErie
#' 
#' @docType data
#' 
#' @format A data frame of 20 observations on the following 2 variables:
#'  \describe{
#'    \item{year}{Year of data.} 
#'    \item{recruits}{Number of recruits (per hour).} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Recruitment time-series 
#'  }
#' 
#' @concept Recruitment
#' @source From Ontario Ministry of Natural Resources, Fish and Wildlife Branch, 1997. Lake Erie fisheries report 1996.  Lake Erie Committee Meeting,
#' Great Lakes Fishery Commission, Ann Arbor, Michigan. pp 26.  Obtained from Ransom Myers online database which was (is?) at http://ram.biology.dal.ca/~myers/data.html.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(RBSmeltErie)
#' str(RBSmeltErie)
#' head(RBSmeltErie)
#' plot(recruits~year,data=RBSmeltErie,type="l")
#' 
NULL