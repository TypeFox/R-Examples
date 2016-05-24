#' @name print.HydeNetwork
#' @export 
#' @method print HydeNetwork
#' 
#' @title Print a HydeNetwork
#' @details Prints a HydeNetwork object with a brief summary of each node.
#' 
#' @param x a \code{HydeNetwork object}
#' @param ... additional arguments to be passed to print methods.  Currently 
#'   none in use.
#'   
#' @details The summary of each node follows the format\cr
#'   node name | parents\cr
#'   node type (parameter)\cr
#'   estimated from data (or not)\cr
#'   formula
#'   
#' @author Jarrod Dalton and Benjamin Nutter
#' @examples
#' data(PE, package="HydeNet")
#' Net <- HydeNetwork(~ wells + 
#'                      pe | wells + 
#'                      d.dimer | pregnant*pe + 
#'                      angio | pe + 
#'                      treat | d.dimer*angio + 
#'                      death | pe*treat) 
#' print(Net)  
#' print(Net, d.dimer) 
#'                     
#' Net <- setNode(Net, d.dimer, 
#'                   nodeType='dnorm', mu=fromData(), tau=fromData(), 
#'                   nodeFormula = d.dimer ~ pregnant + pe,
#'                   nodeFitter='lm')
#' print(Net, d.dimer)
#'     
print.HydeNetwork <- function(x, ...){
  Hyde.nm <- as.character(substitute(x))
   
  #* Requested Nodes
  requested_nodes <- as.character(substitute(list(...)))[-1]
  if (length(requested_nodes) == 0) requested_nodes <- x$nodes
  
  bad_nodes <- requested_nodes[!requested_nodes %in% x$nodes]
  if (length(bad_nodes) > 0)
    stop(paste0("The following nodes are not found in ", substitute(x), ": ", 
                paste(bad_nodes, collapse=", ")))
  
  #* Node Summary Function
  nodeSummary <- function(node){
    nodeName <- if (!is.null(x$parents[[node]]))
                    paste(node, "|", paste(x$parents[[node]], collapse=" * "))
                else node
    
    nodeType <- if (is.null(x$nodeType[[node]])) "Unspecified" else x$nodeType[[node]]
    
    nodeParam <- if (is.null(x$nodeParams[[node]])) "Unspecified" 
                 else{
                   paste(paste(names(x$nodeParams[[node]]), "=", 
                               x$nodeParams[[node]]), collapse=", ")
                 }
    if (nodeType != "Unspecified") nodeType <- paste0(nodeType, "(", nodeParam, ")")
    
    Formula <- paste0(x$nodeFitter[[node]], ": ", deparse(x$nodeFormula[[node]]))
   
    return(paste(nodeName, nodeType, Formula, sep="\n"))
  }
  
  nodeSummaries <- paste(sapply(requested_nodes, 
                                nodeSummary), 
                         collapse="\n\n")

  cat("A Probabilistic Graphical Network", sep=" ")
  cat(paste("\nHas data attached:", if(is.null(x$data)) "No" else "Yes"))
  cat(paste0("\n\n", nodeSummaries))
}
