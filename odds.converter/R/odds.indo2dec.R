#' Convert Indonesian odds to Decimal odds
#'
#' @param x A vector of Indonesian odds
#'
#' @return A vector of Hong Kong odds
#'
#' @examples
#' odds.indo2dec(c(1.93,2.05))
odds.indo2dec <- function (x) {
        ifelse (x > -1 & x < 1,NA,odds.us2dec(odds.indo2us(x)))
}