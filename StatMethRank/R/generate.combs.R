#' Generate all possible combinations of k elements out of n 
#'
#' This function generates all combinations of n elements taken k at a time. It is just 
#' a reference to the combn function in utils package.
#'
#' @param n number of the whole elements
#' @param k number of elements to choose (default 1)
#' @return a matrix or a vector
#' @export
#' @author Scott Chasalow wrote the original in 1994 for S; 
#' R package combinat and documentation by Vince Carey stvjc@@channing.harvard.edu; 
#' small changes by the R core team.
#' @examples
#' generate.combs(10, 6)
#' @seealso utils::combn

generate.combs <- function (n, k = 1)
{
    return(combn(n, k))
}