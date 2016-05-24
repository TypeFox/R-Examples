#' Convert Probabilities to Indonesian odds
#'
#' @param x A vector of Probabilities
#'
#' @return A vector of Indonesian odds
#'
#' @examples
#' odds.prob2indo(c(0.5,0.6))
odds.prob2indo <- function (x){
        ifelse (x <= 0 | x >= 1,NA,odds.us2indo(odds.prob2us(x)))
}