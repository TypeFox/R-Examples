#' Convert Decimal Odds to Probabilities
#'
#' @param x A vector of Decimal odds
#'
#' @return A vector of Probabilities
#'
#' @examples
#' odds.dec2prob(c(1.93,2.05))
odds.dec2prob <- function (x) {
        ifelse(x <=1,NA,1/x)    
}