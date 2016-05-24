#' Convert Indonesian odds to Probabilities
#'
#' @param x A vector of Indonesian odds
#'
#' @return A vector of Probabilities
#'
#' @examples
#' odds.indo2hk(c(1.93,2.05))
odds.indo2prob <- function(x){
        ifelse (x > -1 & x < 1,NA, odds.us2prob(odds.indo2us(x)))
}
