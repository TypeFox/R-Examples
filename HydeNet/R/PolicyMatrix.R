#' @name policyMatrix
#' @export policyMatrix
#' 
#' @title Construct Policy and Decision Matrices
#' @description It may be of interest to compare posterior distributions of 
#'   variables dependent on differing levels of the decision nodes.  For 
#'   example, how might choosing a different engine in a car affect the 
#'   fuel consumption.  \code{policyMatrix} provides a quick utility to 
#'   begin defining the policy matrix on which decisions can be made.
#'   
#' @param network A HydeNetwork object
#' @param ... Named arguments with vectors of the policy values.  
#' @param useDefaultPolicies A logical indicating if the default policy
#'   values from decision nodes should be used for any decision nodes 
#'   not named in \code{...}.  If a node is named in \code{...} and also
#'   has default policy values, the user supplied values take precedence.
#' 
#' @details  When \code{...} is not used, the default policy matrix is defined
#'   as all possible combinations of the levels in the network's decision 
#'   nodes.  If no decision nodes are defined, an error is returned.  Note that
#'   the default policy matrix returns JAGS-ready values, which are numeric 
#'   according to the level number of a factor.  In user-defined matrices, 
#'   character values are supported and will be converted to numerics when the
#'   JAGS model is compiled.
#'   
#'   Semi-custom policy matrices can be defined by providing the values of each
#'   node to be considered.  When manually supplying values, the nodes must 
#'   exist in network, but the requirement that they be decision nodes is not
#'   enforced.  Thus, it is possible to include numeric values in a decision 
#'   matrix.
#'   
#'   Policy matrices can be passed to \code{HydePosterior} to run posterior 
#'   distributions on each row of the policy matrix.  There is nothing 
#'   particularly special about the policy matrices returned by 
#'   \code{policyMatrix}; they are simply data frame that require names drawn
#'   from the nodes in the network.  Any data frame can be passed to 
#'   \code{HydePosterior} and a check is done there to confirm all of the 
#'   column names match a node in the network.
#'   
#'   Whenever a node is identified as a deterministic node, its policy values
#'   are forced to \code{NULL}, regardless of what the user has specified.
#'   
#' @return Returns a data frame built by \code{expand.grid} and intended to be
#'   used with \code{HydePosterior}.
#'   
#' @author Jarrod Dalton and Benjamin Nutter
#' 
#' @examples
#' data(PE, package="HydeNet")
#' Net <- HydeNetwork(~ wells + 
#'                      pe | wells + 
#'                      d.dimer | pregnant*pe + 
#'                      angio | pe + 
#'                      treat | d.dimer*angio + 
#'                      death | pe*treat,
#'                      data = PE) 
#'                      
#' 
#'                  
#' Net <- setDecisionNodes(Net, angio, treat)
#' plot(Net)
#'
#' policyMatrix(Net)
#' 
#' policyMatrix(Net, treat="No", angio = c("No", "Yes"))
#' 
 
policyMatrix <- function(network, ..., useDefaultPolicies = TRUE){
  Check <- ArgumentCheck::newArgCheck()
  utility_nodes <- names(network$nodeUtility)[sapply(network$nodeUtility,
                                                     identity)]
  
  decision_nodes <- names(network$nodeDecision)[sapply(network$nodeDecision, 
                                                       identity)]
  decision_nodes <- decision_nodes[!decision_nodes %in% utility_nodes]
  
  defaultPolicies <- network$nodePolicyValues[decision_nodes]
  
  policies <- list(...)
  policies <- policies[!names(policies) %in% utility_nodes]
  
  if (length(policies) < 1) return(defaultPolicyMatrix(network))
  
  if (!(all(names(policies) %in% network$nodes))){
    not_in_network <- names(policies)[!names(policies) %in% network$nodes]
    ArgumentCheck::addError(paste0("The following input nodes do not exist in '",
                                   substitute(network), "': ",
                                   paste(not_in_network, collapse=", "), "."),
                            Check)
  }
  
  ArgumentCheck::finishArgCheck(Check)
  
  if (any(names(policies) %in% names(defaultPolicies))){
    defaultPolicies <- defaultPolicies[!names(defaultPolicies) %in% names(policies)]
  }
  
  expand.grid(c(policies, defaultPolicies), stringsAsFactors=FALSE)
}

#' @rdname policyMatrix

defaultPolicyMatrix <- function(network){
  Check <- ArgumentCheck::newArgCheck()
  utility_nodes <- names(network$nodeUtility)[sapply(network$nodeUtility,
                                                     identity)]
  
  decision_nodes <- names(network$nodeDecision)[sapply(network$nodeDecision, 
                                                       identity)]
  
  decision_nodes <- decision_nodes[!decision_nodes %in% utility_nodes]
  
  if (length(decision_nodes) == 0)
  ArgumentCheck::addError(paste0("There are no decision nodes in '", substitute(network), "'."),
                          Check)
  
  decision_options <- lapply(decision_nodes, policyMatrixValues, network)
  names(decision_options) <- decision_nodes
  
  expand.grid(decision_options) 
}
