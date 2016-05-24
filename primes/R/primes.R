#' @title Functions for Identifying and Generating Prime Numbers
#' @name primes
#' @description blardfdsfsdfsd
#'
#' @useDynLib primes
#' @importFrom Rcpp sourceCpp
#' @docType package
#' @aliases primes primes-package
NULL

#' @title Generate and Test for Prime Numbers
#' @description generate prime numbers or test whether a sequence of numbers you have are
#' prime or not.
#'
#' @param x an integer vector containing elements you want to determine the primality of.
#'
#' @param min the value to generate primes from.
#'
#' @param max the maximum value to generate prime numbers up to.
#'
#' @details \code{is_prime} and \code{generate_primes} rely on Wilson's theorem to test for a number's primality;
#' as primality algorithms go, this is actually a very \emph{slow} approach - in theory. In practice, because of the
#' limits R institutes around integer sizes, it's fast enough for our needs. For example, 10m numbers, all 2^30-sized,
#' can be tested for primality using this package in 100ms.
#'
#'
#' @examples
#'
#' #Test for primality
#' is_prime(1299827)
#' # [1] TRUE
#'
#' generate_primes(max =12)
#' # [1]  2  3  5  7 11
#' @rdname prime
#' @name prime
#' @export
is_prime <- function(x){
  is_prime_vector(x)
}

#'@rdname prime
#'@export
generate_primes <- function(min = 0, max){
  generate_primes_(min, max)
}
