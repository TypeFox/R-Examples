#' @title Utilities for distributions supported in this package
#' @name distname-utils
#' @aliases check_distname get_distnames get_distname_family
#' 
#' @description
#' The Lambert W\eqn{\times} F framework can take any (continuous) random variable with distribution
#' F and make it skewed (\code{type = "s"}), heavy tailed (\code{type = "h"}),
#' or both (\code{type = "hh"}).
#' 
#' In principle, this works for any F.  Of course, this package implements only a finite
#' number of distributions, which can be specified with the \code{distname} argument.
#' Most functions in this package, however, also allow you to pass your own distribution and parameters
#' and create a Lambert W\eqn{\times} F version of it.
#' 
#' @seealso \code{\link{create_LambertW_input}}, \code{\link{create_LambertW_output}}.
#' 
#' @inheritParams common-arguments
#' @keywords misc
#' 
NULL