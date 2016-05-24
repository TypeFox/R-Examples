#' datafsm: A package for estimating FSM models.
#'
#'It relies on the \strong{GA} package:
#'Luca Scrucca (2013). GA: A Package for Genetic Algorithms in R. 
#'Journal of Statistical Software, 53 (4), 1-37. 
#'URL http://www.jstatsoft.org/v53/i04/.
#'
#' @section datafsm functions:
#' \code{datafsm}'s main function for estimating a fsm decision
#' model:
#' \enumerate{
#' \item \code{\link{evolve_model}}
#' }
#' \code{datafsm}'s helper functions:
#' \enumerate{
#' \item \code{\link{evolve_model_cv}}
#' \item \code{\link{var_imp}}
#' \item \code{\link{decode_state_mat}}
#' \item \code{\link{decode_action_vec}}
#' \item \code{\link{fitnessCPP}}
#' \item \code{\link{build_bitstring}}
#' \item \code{\link{compare_fsm}}
#' }
#' 
#' @docType package
#' @name datafsm
NULL

#' @importFrom methods setClass setGeneric setMethod setRefClass
NULL
