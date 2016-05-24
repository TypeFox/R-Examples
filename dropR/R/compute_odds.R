#' Compute Odds From Probabilities
#' 
#' Compute odds from probabilities. The function is vectorized and can handle a vector of probabilities. 
#' 
#' @param p vector of probabilities. May not be larger than 1 or smaller than zero.
#' @examples
#' get_odds(.8) # 4
#' @export
get_odds <- function(p){
  if(!all(p <= 1 & p >= 0)) stop('Input is not a probability!')
  p / (1-p)
} 

