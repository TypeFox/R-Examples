#' Convert Indonesian odds to Malaysian odds
#'
#' @param x A vector of Indonesian odds
#'
#' @return A vector of Malaysian odds
#'
#' @examples
#' odds.indo2malay(c(1.93,2.05))
odds.indo2malay <- function (x) {
        ifelse (x > -1 & x < 1,NA,odds.us2malay(odds.indo2us(x)))
}