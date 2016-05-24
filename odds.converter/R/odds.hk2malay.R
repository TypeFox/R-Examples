#' Convert Hong Kong Odds to Malaysian Odds
#'
#' @param x A vector of Hong Kong odds
#'
#' @return A vector of Malaysian odds
#'
#' @examples
#' odds.hk2malay(c(1.93,0.05))
odds.hk2malay <- function (x){
        ifelse (x <= 0,NA,odds.us2malay(odds.hk2us(x)))
}