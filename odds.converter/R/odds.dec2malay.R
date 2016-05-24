#' Convert Decimal Odds to Malaysian odds
#'
#' @param x A vector of Decimal odds
#'
#' @return A vector of Malaysian odds
#'
#' @examples
#' odds.dec2malay(c(1.93,2.05))
odds.dec2malay <- function (x) {
        ifelse (x <= 1,NA,odds.us2malay(odds.dec2us(x)))
}
