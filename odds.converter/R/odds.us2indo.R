#' Convert US odds to Indonesian odds
#'
#' @param x A vector of US odds
#'
#' @return A vector of Indonesian odds
#'
#' @examples
#' odds.us2indo(c(-200,150))
odds.us2indo <- function (x) {
        ifelse(x > -100 & x <100,NA,x/100)   
}