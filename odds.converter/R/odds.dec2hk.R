#' Convert Decimal Odds to Hong Kong odds
#'
#' @param x A vector of Decimal odds
#'
#' @return A vector of Hong Kong odds
#'
#' @examples
#' odds.dec2hk(c(1.93,2.05))
odds.dec2hk <- function (x) {
        ifelse (x <= 1,NA,x-1)
}
