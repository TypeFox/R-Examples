#' Convert Indonesian odds to Hong Kong odds
#'
#' @param x A vector of Indonesian odds
#'
#' @return A vector of Hong Kong odds
#'
#' @examples
#' odds.indo2hk(c(1.93,2.05))
odds.indo2hk <- function (x) {
        ifelse (x > -1 & x < 1,NA,odds.us2hk(odds.indo2us(x)))
}