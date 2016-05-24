#' Convert Malaysian odds to Hong Kong odds
#'
#' @param x A vector of Malaysian odds
#'
#' @return A vector of Hong Kong odds
#'
#' @examples
#' odds.malay2hk(c(1.93,2.05))
odds.malay2hk <- function (x) {
        ifelse (x < -1 | x > 1,NA,odds.us2hk(odds.malay2us(x)))
}