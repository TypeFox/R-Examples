#' Convert Hong Kong Odds to US Odds
#'
#' @param x A vector of Hong Kong odds
#'
#' @return A vector of US odds
#'
#' @examples
#' odds.hk2us(c(1.93,0.05))
odds.hk2us <- function (x){
        ifelse (x <= 0,NA,ifelse(x > 1, 100*x,-100/x))
}