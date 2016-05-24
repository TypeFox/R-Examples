#' @title Recruitment time-series for Walleye in Lake Erie, 1959-1972.
#' 
#' @description Walleye (\emph{Sander vitreus}) recruitment time-series for Lake Erie, 1959-1972.
#' 
#' @name WalleyeErie
#' @docType data
#' 
#' @format A data frame of 14 observations on the following 6 variables:
#' \describe{
#'  \item{year}{Year of data.} 
#'  \item{recruits}{Number of recruits (per 1000 ft of net).} 
#' }
#' 
#' @section Topic(s): \itemize{
#'  \item Recruitment time-series 
#' }
#' 
#' @concept Recruitment
#' 
#' @source From Wolfert, D. R. 1981. The commercial fishery for walleyes in New York waters of Lake Erie, 1959-1978. North American Journal of Fisheries Management 1: 112-126.  Obtained from Ransom Myers online database which was (is?) at http://ram.biology.dal.ca/~myers/data.html.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WalleyeErie)
#' str(WalleyeErie)
#' head(WalleyeErie)
#' plot(recruits~year,data=WalleyeErie,type="l")
#' 
NULL