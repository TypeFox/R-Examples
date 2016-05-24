#' Convert Malaysian odds to Indonesian odds
#'
#' @param x A vector of Malaysian odds
#'
#' @return A vector of Indonesian odds
#'
#' @examples
#' odds.malay2indo(c(1.93,2.05))
odds.malay2indo <- function (x) {
        ifelse (x < -1 | x > 1,NA,odds.us2indo(odds.malay2us(x)))
}