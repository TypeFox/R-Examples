#' @name Domains
#' @aliases Domains
#' @title Domains
#' 
#' @description Given a main vector, an auxiliary one and a value of the latter, identifies positions of the auxiliary vector corresponding
#'  to values other than the given one. Then, turns zero values of the main vector corresponding to these positions.
#' 
#' @usage Domains (y, domains, value)
#' @param y A numeric main vector of size n
#' @param domains A numeric/character/logic auxiliary vector of size n
#' @param value A value of the auxiliary vector
#' @return A numeric vector, copy of \code{y}, with some values turned zero depending on values of \code{domains} and \code{value}
#' @examples
#' ##########   Example 1   ##########
#' U <- c(13, 18, 20, 14, 9)
#' #Let build an auxiliary vector indicating whether values in U are above or below the mean.
#' aux <- c("Below", "Above", "Above", "Below", "Below")
#' #Now, only values below the mean remain, the other ones are turned zero.
#' Domains (U, aux, "Below")
#' 
#' ##########   Example 2   ##########
#' data(DatA)
#' attach(DatA)
#' #Let calculate total feeding expenses corresponding to households in domain a.
#' sum (Domains (Feed, Domain, "a"))
#' @export
Domains = function (y, domains, value)
{
	if (any(is.na(y))) 
        	stop("There are missing values in y.")
    	if (length(y) != length(domains)) 
        	stop("y and domains vector have different sizes.")
	
	return (y * (domains == value))
}