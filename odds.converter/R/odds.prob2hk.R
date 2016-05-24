#' Convert Probabilities to Hong Kong odds
#'
#' @param x A vector of Probabilities
#'
#' @return A vector of Hong Kong odds
#'
#' @examples
#' odds.prob2hk(c(0.5,0.6))
odds.prob2hk <- function (x){
        ifelse (x <= 0 | x >= 1,NA,(1/x)-1)
}