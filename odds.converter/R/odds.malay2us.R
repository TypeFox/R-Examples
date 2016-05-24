#' Convert Malaysian odds to US odds
#'
#' @param x A vector of Malaysian odds
#'
#' @return A vector of US odds
#'
#' @examples
#' odds.malay2us(c(0.5,-0.6))
odds.malay2us <- function (x) {
        ifelse (x < -1 | x > 1,NA,-100/x)
}