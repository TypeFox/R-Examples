#' @name setPolicyValues
#' @export
#' 
#' @title Assign Default Policy Values
#' @description By default, \code{HydeNet} uses factor levels for policy
#'   values in a decision node, assuming the decision node is a factor 
#'   variable. In cases where the decision node is a numeric variable, 
#'   \code{HydeNet} will first try to assign the first, second, and 
#'   third quartiles as policy values.  \code{setPolicyValues} allows 
#'   the user flexibility in which values are actually used in the 
#'   decision network.  It can also be used to restrict the levels of 
#'   a factor variable to a subset of all levels.  Policy values may also 
#'   be set in \code{setNode}, but \code{setPolicyValues} makes it possible
#'   to set the values for multiple nodes in one call.
#'   
#' @param network A Hyde Network object
#' @param ... arguments named for nodes in the network.  The value of each
#'   argument will be assigned to the \code{nodePolicyValues} element of 
#'   the \code{HydeNetwork} object.
#'   
#' @author Jarrod Dalton and Benjamin Nutter

setPolicyValues <- function(network, ...){
  policyValues <- list(...)
  for (i in names(policyValues)){
    network$nodePolicyValues[[i]] <- policyValues[[i]]
  }
  network
}