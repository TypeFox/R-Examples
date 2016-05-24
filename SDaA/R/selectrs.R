#' Steps used in Selecting an SRS  
#'
#' Steps used in selecting the simple random sample (SRS) in Example 2.4
#' of Lohr(1999).
#' @name selectrs
#' @docType data
#' @format Data frame with the following 5 variables: 
#' \describe{
#'   \item{a}{random number generated between 0 and 1}
#'   \item{b}{ceiling(3048*RN), with RN the random number
#'            in column \code{a}}
#'   \item{c}{distinct values in column \code{b}}
#'   \item{d}{new values generated to replace duplicates in \code{b}}
#'   \item{e}{final set of distinct values to be used in sample}
#' }
#' @note the set of indices in column \code{e} was used to select
#' observations from \code{agpop} into dataset \code{agsrs}.
#' @references Lohr (1999). Sampling: Design and Analysis, Duxbury, p. 31--34 
#' and 444. 
#' @export
NULL
