#' Convert US odds to Malaysian odds
#'
#' @param x A vector of US odds
#'
#' @return A vector of Malaysian odds
#'
#' @examples
#' odds.us2malay(c(-200,150))
odds.us2malay <- function (x) {
        ifelse(x > -100 & x < 100,NA,-100/x)
}