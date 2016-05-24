#' Convert Indonesian odds to US odds
#'
#' @param x A vector of Indonesian odds
#'
#' @return A vector of US odds
#'
#' @examples
#' odds.indo2us(c(1.93,2.05))
odds.indo2us <- function (x) {
        ifelse (x > -1 & x < 1,NA,100 * x)
}