#' Convert Probabilities to Decimal odds
#'
#' @param x A vector of Probabilities
#'
#' @return A vector of Decimal odds
#'
#' @examples
#' odds.prob2dec(c(0.5,0.6))
odds.prob2dec <- function (x) {
        ifelse (x <= 0 | x >= 1,NA,1/x)
}
                