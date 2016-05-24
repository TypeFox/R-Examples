#' Convert Probabilities to US odds
#'
#' @param x A vector of Probabilities
#'
#' @return A vector of US odds
#'
#' @examples
#' odds.prob2us(c(0.5,0.6))
odds.prob2us <- function (x) {
        ifelse (x <= 0 | x >= 1,NA,ifelse(x <=0.5,100 * (1/x - 1), -100/(1/x - 1)))
}