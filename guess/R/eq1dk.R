#' Constraints: Sum to 1
#' @description Constraints that some params sum to 1. Used Internally. For data with DK.
#' Functions for constraining lambdas to sum to 1 and to bound params between 0 and 1
#' @keywords internal
#' 
#' @param x   lgg, lgk, lgc, lkk, lcg, lck, and lcc
#' @param g1  guess
#' @param data transition matrix

eqn1dk = function(x, g1=NA, data) {
	sum(x[1:7])
}
