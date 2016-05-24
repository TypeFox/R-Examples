#' Convert Decimal Odds to all other formats
#'
#' @param x A vector of Decimal odds
#'
#' @return A dataframe 
#'
#' @examples
#' odds.dec2all(c(1.93,2.05))
odds.dec2all <- function(x) {
  data.frame(  DECIMAL = x,
               US =  round(odds.dec2us(x),4),
               PROB = round(odds.dec2prob(x),4),
               HK = round(odds.dec2hk(x),4),
               INDO = round(odds.dec2indo(x),4),
               MALAY = round(odds.dec2malay(x),4)
               )
}

#' Convert US Odds to all other formats
#'
#' @param x A vector of US odds
#'
#' @return A dataframe 
#'
#' @examples
#' odds.us2all(c(-200,105))
odds.us2all <- function(x) {
  data.frame(  US = x,
               DECIMAL =  round(odds.us2dec(x),4),
               PROB = round(odds.us2prob(x),4),
               HK = round(odds.us2hk(x),4),
               INDO = round(odds.us2indo(x),4),
               MALAY = round(odds.us2malay(x),4)
               )
}

#' Convert Probabilities to  to all other formats
#'
#' @param x A vector of Probabilities
#'
#' @return A dataframe
#'
#' @examples
#' odds.prob2all(c(0.5,0.6))
odds.prob2all <- function(x) {
  data.frame(  PROB = x,
               DECIMAL =  round(odds.prob2dec(x),4),
               US = round(odds.prob2us(x),4),
               HK = round(odds.prob2hk(x),4),
               INDO = round(odds.prob2indo(x),4),
               MALAY = round(odds.prob2malay(x),4)
               )
}


#' Convert Hong Kong Odds to all other formats
#'
#' @param x A vector of Hong Kong odds
#'
#' @return A dataframe
#'
#' @examples
#' odds.hk2all(c(1.93,0.05))
odds.hk2all <- function(x) {
  data.frame(  HK = x,
               DECIMAL =  round(odds.hk2dec(x),4),
               US = round(odds.hk2us(x),4),
               PROB = round(odds.hk2prob(x),4),
               INDO = round(odds.hk2indo(x),4),
               MALAY = round(odds.hk2malay(x),4)
               )
}


#' Convert Indonesian to all other formats
#'
#' @param x A vector of Indonesian odds
#'
#' @return A vector of US odds
#'
#' @examples
#' odds.indo2all(c(1.93,2.05))
odds.indo2all <- function(x) {
  data.frame(  INDO = x,
               DECIMAL =  round(odds.indo2dec(x),4),
               US = round(odds.indo2us(x),4),
               PROB = round(odds.indo2prob(x),4),
               HK = round(odds.indo2hk(x),4),
               MALAY = round(odds.indo2malay(x),4)
               )
}

#' Convert Malaysian Odds to all other formats
#'
#' @param x A vector of Malaysian odds
#'
#' @return A dataframe
#'
#' @examples
#' odds.malay2all(c(0.5,-0.6))
odds.malay2all <- function(x) {
  data.frame(  MALAY = x,
               DECIMAL =  round(odds.malay2dec(x),4),
               US = round(odds.malay2us(x),4),
               PROB = round(odds.malay2prob(x),4),
               HK = round(odds.malay2hk(x),4),
               INDO = round(odds.malay2indo(x),4)
  )
}
