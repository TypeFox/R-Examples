#' @title Utilities for transformation vector tau
#' @name tau-utils
#' @aliases check_tau tau2theta complete_tau get_initial_tau tau2type
#' 
#' @description All functions here are for the transformation parameter vector
#'     \eqn{\tau = (\mu_x, \sigma_x, \gamma, \delta, \alpha)}.
#' 
#' @inheritParams common-arguments
#' @return
#' \code{check_tau} throws an error if \eqn{\tau} does not define a proper
#'     transformation.
#' 
#' \code{complete_tau} returns a named numeric vector.
#' 
#' \code{get_initial_tau} returns a named numeric vector.
#' 
#' \code{tau2theta} returns a list with entries \code{alpha}, \code{beta},
#'     \code{gamma}, and \code{delta}.
#' 
#' \code{tau2type} returns a string: either \code{"s"}, \code{"h"}, or
#'     \code{"hh"}.
#' 
#' @keywords math utilities
#' 
#' 
NULL
