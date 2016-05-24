#' Convert Malaysian odds to Probabilities
#'
#' @param x A vector of Malaysian odds
#'
#' @return A vector of Probabilities
#'
#' @examples
#' odds.malay2prob(c(1.93,2.05))
odds.malay2prob <- function (x) {
        ifelse(x < -1 | x > 1,NA,odds.us2prob(odds.malay2us(x)))
}