#' Convert Hong Kong Odds to Indonesian Odds
#'
#' @param x A vector of Hong Kong odds
#'
#' @return A vector of Indonesian odds
#'
#' @examples
#' odds.hk2indo(c(1.93,0.05))
odds.hk2indo <- function (x){
        ifelse (x <= 0,NA,odds.us2indo(odds.hk2us(x)))
}