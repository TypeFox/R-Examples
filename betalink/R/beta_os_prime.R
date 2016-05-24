#' @title Measure the distance between a network and its metaweb
#' @description
#' Returns the values of beta OS', i.e. the distace between all realizations, and the revelant subset from the metaweb
#'
#' @param N a list of networks
#' @param ... additional arguments to be passed to \link{betalink}
#'
#' @return An array of the values of Beta OS'
#'
#' @export
beta_os_prime <- function(N, ...){
	M <- metaweb(N)
	os_prime <- plyr::laply(N, function(x) betalink(x, M, ...)$OS)
   names(os_prime) <- names(N)
	return(os_prime)
}
