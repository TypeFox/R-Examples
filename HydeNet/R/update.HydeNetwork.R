#' @name update.HydeNetwork
#' @export 
#' @method update HydeNetwork
#' 
#' @title Update Probabilistic Graphical Network
#' @description Add or remove nodes or add parents within a \code{HydeNetwork}
#'   model.
#'   
#' @param object A \code{HydeNetwork} object
#' @param formula A formula statement indicating the changes to the network.
#' @param ... Additional arguments to be passed to other methods.  Current,
#'   none are used.
#'   
#' @details Adding or removing nodes is fairly straightforward if you are 
#'   removing a complete node (along with its parents).  If you just wish to 
#'   remove a parent, this doesn't work well yet.  I will have to work out 
#'   a solution to remove parent relationships individually.  I had hoped
#'   there would be an \code{update.dag} method, but no luck.  This will continue
#'   to be developed in the future, but the arguments will remain the same.
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
#'                      death | pe*treat)
#'                      
#' plot(Net)
#' 
#' Net <- update(Net, . ~ . - pregnant)
#' plot(Net)
#'  
                   
update.HydeNetwork <- function(object, formula, ...){
  new_formula <- rewriteHydeFormula(object$network_formula, formula)
  
  NEW <- HydeNetwork(new_formula, data=object$data)
  
  lostParents <- lapply(names(NEW$parents),
         function(nm){
           setdiff(object$parents[[nm]], NEW$parents[[nm]])
         })
  names(lostParents) <- names(NEW$parents)
  
  if (any(sapply(lostParents, length) > 0)){
    lostParents <- lostParents[sapply(lostParents, length) > 0]
    warning(paste0("The following nodes lost parents in the update--please redefine the node formula:\n",
                   paste0("    ", names(lostParents), ": ", sapply(lostParents, paste, collapse=", "),
                          collapse="\n")))
  }
  
  
  
  
  NEW$nodeType[names(object$nodeType)] <- object$nodeType
  NEW$nodeFormula[names(object$nodeFormula)] <- object$nodeFormula
  NEW$nodeFitter[names(object$nodeFitter)] <- object$nodeFitter
  NEW$nodeFitterArgs[names(object$nodeFitterArgs)] <- object$nodeFitterArgs
  NEW$nodeParams[names(object$nodeParams)] <- object$nodeParams
  NEW$nodeData[names(object$nodeData)] <- object$nodeData
  
  return(NEW)  
}
