#' Convert US odds to Decimal odds
#'
#' @param x A vector of US odds
#'
#' @return A vector of Decimal odds
#'
#' @examples
#' odds.us2dec(c(-200,150))
odds.us2dec <- function (x) {
        ifelse (x >= 100 & x < 100,NA,ifelse(x >0,x/100+1,-100/x+1))
}              
