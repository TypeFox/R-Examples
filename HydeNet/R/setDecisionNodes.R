#' @name setDecisionNodes
#' @export
#' 
#' @title Classify Multiple Nodes as Decision or Utility Nodes
#' @description Depending on how your Hyde Network was built, you may not have
#'   had the opportunity to declare which nodes are decision or utlity nodes.  
#'   For instance,
#'   when passing a list of models, there is no way to indicate in the model
#'   object that the node should be considered a decision node.  As a matter 
#'   of convenience, these function will set any nodes indicated to decision 
#'   or utility nodes.  It will make no other modifications to a node's definition.
#'   
#' @param network A Hyde Network object
#' @param ... Networks to be classified as decision nodes.  These may be quoted
#'   or unquoted.
#'   
#' @author Jarrod Dalton and Benjamin Nutter
#' 

setDecisionNodes <- function(network, ...){
  nodes <- as.character(substitute(list(...)))[-1]
  network$nodeDecision[nodes] <- lapply(network$nodeDecision[nodes], function(x) TRUE)
  network
}

#' @rdname setDecisionNodes
#' @export

setUtilityNodes <- function(network, ...){
  nodes <- as.character(substitute(list(...)))[-1]
  network$nodeUtility[nodes] <- lapply(network$nodeUtility[nodes], function(x) TRUE)
  network
}


