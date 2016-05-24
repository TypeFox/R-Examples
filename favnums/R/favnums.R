#' @title A Dataset of Favourite Numbers
#' @name favnums
#' @description This package provides a dataset of favourite numbers, selected from an online poll of over
#' 30,000 people by Alex Bellos
#' 
#' @seealso \code{\link{favourite_numbers}} and the \href{http://pages.bloomsbury.com/favouritenumber}{original dataset}.
#' @docType package
#' @aliases favnums favnums-package
NULL

#' @title Favourite Numbers based on an online poll
#'
#' @description A dataset containing the favourite numbers selected by over
#' 30,000 people in an online poll.
#'
#' @format A data frame with 1123 rows and 4 variables:
#' \describe{
#'   \item{number}{the actual number chosen. May be NA in the case of "imaginary" numbers, or Infinite.}
#'   \item{frequency}{the number of times this number was chosen.}
#'   \item{percentage}{the percentage of user answers a particular number represents.}
#'   \item{description}{descriptions of the number's importance, as provided by Alex Bellos. Often NA.}
#' }
#' @source \url{http://pages.bloomsbury.com/favouritenumber}
#' 
#' @examples
#' head(favourite_numbers)
"favourite_numbers"
