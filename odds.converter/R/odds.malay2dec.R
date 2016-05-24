#' Convert Malaysian odds to Decimal odds
#'
#' @param x A vector of Malaysian odds
#'
#' @return A vector of Decimal odds
#'
#' @examples
#' odds.malay2dec(c(0.5,-0.6))
odds.malay2dec <- function (x) {
        ifelse (x < -1 | x > 1,NA,odds.us2dec(odds.malay2us(x)))
}