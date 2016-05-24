#' Calculate multinomial coefficient
#'
#' This function calculates the multinomial coefficient
#' \deqn{\frac{(\sum n_j)!}{\prod n_j!}.}{(\sum n_j)! / \prod n_j!.}
#' where \eqn{n_j}{n_j}'s are the number of multiplicities in the multiset.
#'
#' The result is not reliable if the retunred value is too large (around 2^53) due to limition of integers.
#' @param n a vector of group sizes
#' @return multinomial coefficient
#' @examples
#' # (3+1+1)!/ (3! 1! 1!) = 20
#' multichoose(c(3,1,1))
#' @export
multichoose <- function(n){
    round(exp(lgamma(sum(n) + 1) - sum(lgamma(n + 1))))
}


#' Calculate the number of r-permutations of a multiset
#' @param f the frequencies of the mutliset
#' @param r the number of object drawn from the multiset
#' @return the number of r-permutations
#' @examples
#' x = c("a","a","b")
#' # possible permutations of size 2 are "aa", "ab" and "ba".
#' np_multiset(table(x), 2) # = 3
#' @export
#' @import polynom
np_multiset <- function(f, r){
    p <- polynomial(1)
    alpha <- lgamma(r + 1) / length(f)
    for (i in f){
        j <- pmin(i, r)
        p <- p * polynomial(exp(c(alpha, alpha - lgamma( (1:j) + 1 ))))
        p <- polynomial(p[1:min(length(p), r + 1)])
    }
    return(round(p[r + 1]))
}


#' Calculate the number of r-combinations of a multiset
#' @param f the frequencies of the mutliset
#' @param r the number of object drawn from the multiset
#' @return the number of combinations
#' @examples
#' x <- c("a","a","b")
#' # possible combinations of size 2 are "aa" and "ab".
#' nc_multiset(table(x), 2) # <- 2
#' @export
#' @import polynom
nc_multiset <- function(f, r){
    p <- polynomial(1)
    for (i in f){
        p <- p * polynomial(rep.int(1, i + 1))
        p <- polynomial(p[1:min(length(p), r + 1)])
    }
    return(p[r + 1])
}

#' Wrap iterpc objects by iterators::iter
#' @param I the iterpc object
#' @param d number of permutation(s)/combination(s) wanted in each iteration, default to 1
#' @return a iter object compatible with iterators package
#' @import iterators
#' @examples
#' library(iterators)
#' I <- iterpc(5, 2)
#' it <- iter_wrapper(I)
#' nextElem(it)
#' nextElem(it)
#'
#' library(foreach)
#' I <- iterpc(5, 2)
#' it <- iter_wrapper(I)
#' foreach(x=it, .combine=c) %do% { sum(x) }
#' @export
iter_wrapper <- function(I, d=1){
    iter(function() {
        out <- getnext(I, d)
        !is.null(out) || stop("StopIteration", call. = FALSE)
        out
    })
}
