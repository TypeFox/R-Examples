#' Convert Probabilities to Malaysian odds
#'
#' @param x A vector of Probabilities
#'
#' @return A vector of Malaysian odds
#'
#' @examples
#' odds.prob2malay(c(0.5,0.6))
odds.prob2malay <- function (x){
        ifelse (x <= 0 | x >= 1,NA,odds.us2malay(odds.prob2us(x)))
}
