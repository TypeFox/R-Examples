#' Convert US odds to Probabilities
#'
#' @param x A vector of US odds
#'
#' @return A vector of Probabilities
#'
#' @examples
#' odds.us2prob(c(-200,150))
odds.us2prob <- function (x) {
        ifelse (x > -100 & x < 100,NA,ifelse(x >0,100/(100 + x), x/(x - 100)))
}