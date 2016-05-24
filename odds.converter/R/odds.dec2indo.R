#' Convert Decimal Odds to Indonesian odds
#'
#' @param x A vector of Decimal odds
#'
#' @return A vector of Indonesian odds
#'
#' @examples
#' odds.dec2indo(c(1.93,2.05))
odds.dec2indo <- function (x) {
        ifelse (x <= 1,NA,odds.us2indo(odds.dec2us(x)))
}
