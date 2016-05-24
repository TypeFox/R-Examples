#' Convert Hong Kong Odds to Probabilities
#'
#' @param x A vector of Hong Kong odds
#'
#' @return A vector of Probabilities
#'
#' @examples
#' odds.hk2us(c(1.93,0.05))
odds.hk2prob<- function (x){
        ifelse (x <= 0,NA,odds.dec2prob(x+1))
}