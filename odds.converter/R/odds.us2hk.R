#' Convert US odds to Hong Kong odds
#'
#' @param x A vector of US odds
#'
#' @return A vector of Hong Kong odds
#'
#' @examples
#' odds.us2hk(c(-200,150))
odds.us2hk <- function (x){
        odds.us2dec(x)-1
}