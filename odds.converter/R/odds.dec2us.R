#' Convert Decimal Odds to US Odds
#'
#' @param x A vector of Decimal odds
#'
#' @return A vector of US odds
#'
#' @examples
#' odds.dec2us(c(1.93,2.05))
odds.dec2us <- function (x) {
        ifelse (x <= 1,NA,ifelse(x >2, 100 * (x - 1),-100/(x - 1)))
}
