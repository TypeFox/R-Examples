#' @name writeJagsModel
#' 
#' @title Write a Node's JAGS Model
#' @description Constructs the JAGS code that designates the model for the 
#'   node conditioned on its parents.  The parameters for the model may 
#'   be user supplied or estimated from a given data set.
#'  
#' @param network A network of class HydeNetwork
#' @param node A node within \code{network}
#' 
#' @details The manipulations are performed on the \code{nodeParams} element
#'   of the \code{Hyde} network.  A string of JAGS code is returned suitable
#'   for inclusion in the Bayesian analysis.  
#'   
#'   The function will (eventually) travel through a serious of \code{if} 
#'   statements until it finds the right node type.  It will then match
#'   the appropriate arguments to the inputs based on user supplied values or
#'   estimating them from the data.
#'   
#' @author Jarrod Dalton and Benjamin Nutter
#' @seealso \code{\link{writeJagsFormula}}
#' 
#' @examples
#' #* NOTE: writeJagsModel isn't an exported function
#' data(PE, package='HydeNet')
#' Net <- HydeNetwork(~ wells + 
#'                      pe | wells + 
#'                      d.dimer | pregnant*pe + 
#'                      angio | pe + 
#'                      treat | d.dimer*angio + 
#'                      death | pe*treat,
#'                      data = PE)
#' HydeNet:::writeJagsModel(Net, 'pe')
#' HydeNet:::writeJagsModel(Net, 'treat')
#' 

writeJagsModel <- function(network, node){
  node_str <- if (is.character(node)) node else as.character(substitute(node))

  node_params <- network$nodeParams[[node_str]]

  params <- eval(substitute(expectedParameters(network, node, TRUE)))
  
  #*** Type 'determ' (Deterministic Nodes)
  if (network$nodeType[[node_str]] == "determ"){
    if (fromFormula() %in% node_params){
      define <- as.character(network$nodeFormula[[node_str]])
      define <- paste0(define[2], " <- ", define[3])
      model_code <- define
    }
    else model_code <- paste0(node_str, " <- ", network$nodeParams[[node_str]]$define)
  }
  
  #*** Type 'dbern'
  else if (network$nodeType[[node_str]] == "dbern"){
    if (any(c(fromData()) %in% node_params))
      fit <- do.call(network$nodeFitter[[node_str]],
                     c(list(formula = network$nodeFormula[[node_str]],
                            data = if (is.null(network$nodeData[[node_str]])) network$data else network$nodeData[[node_str]]),
                       network$nodeFitterArgs[[node_str]]))
    
    if (node_params['p'] %in% c(fromData(), fromFormula())){
      if (node_params['p'] == fromData())
        node_params['p'] <- writeJagsFormula(fit, network$nodes)
      else if (node_params['p'] == fromFormula())
        node_params['p'] <- rToJags(network$nodeFormula[[node_str]])
    
      node_params['p'] <- as.character(as.formula(node_params[['p']]))[-(1:2)]
    }
    
    model_code <- paste0(node_str, " ~ ",
                         network$nodeType[[node_str]],
                         "(",
                         paste(node_params, collapse=", "),
                         ")")

  }
  
  #*** Type 'dcat'
  else if (network$nodeType[[node_str]] == "dcat"){
    if (!is.null(network$nodeFitter[[node_str]]) && network$nodeFitter[[node_str]] == "cpt"){
      parents <- network$parents[[node_str]]
      bern_parent <- sapply(parents, function(p) network$nodeType[[p]] == "dbern")
      parents[bern_parent] <- paste0("(", parents[bern_parent], "+1)")
      
      model_code <- paste0("pi.", node_str, " <- cpt.", node_str, "[",
                           paste0(parents, collapse=", "), ", ]\n",
                           "   ", node_str, " ~ dcat(pi.", node_str, ")")

    } else if (fromData() %in% node_params){
      node_params["pi"] <- paste0("pi.", node_str)
      pi <- do.call(network$nodeFitter[[node_str]],
                    list(formula = network$nodeFormula[[node_str]],
                         data=if (is.null(network$nodeData[[node_str]])) network$data else network$nodeData[[node_str]]))
      pi <- writeJagsFormula(pi, network$nodes)
      model_code <- 
        c(pi,  
          paste0(node_str, " ~ ",
                 network$nodeType[[node_str]],
                 "(",
                 paste(node_params, collapse=", "),
                 ")"))
      }
      else{
        model_code <- paste0(network$nodeParams[[node_str]]['pi'], "\n   ",
                              node_str, " ~ ", network$nodeType[[node_str]],
                              "(pi.", node_str, ")")
    }
  }
  
  #*** Type 'dnorm'
  else if (network$nodeType[[node_str]] == "dnorm"){
    if (any(c(fromData()) %in% node_params))
      fit <- do.call(network$nodeFitter[[node_str]],
                     c(list(formula = network$nodeFormula[[node_str]],
                          data = if (is.null(network$nodeData[[node_str]])) network$data else network$nodeData[[node_str]]),
                       network$nodeFitterArgs[[node_str]]))
    
    if (node_params['mu'] %in% c(fromData(), fromFormula())){
      if (node_params['mu'] == fromData())
        node_params['mu'] <- writeJagsFormula(fit, network$nodes)
      
      else if (node_params['mu'] == fromFormula())
        node_params['mu'] <- rToJags(network$nodeFormula[[node_str]])
    
      
      node_params['mu'] <- as.character(as.formula(node_params[['mu']]))[-(1:2)]

    }
    
    if (node_params['tau'] %in% c(fromData(), fromFormula())){
      if (node_params['tau'] == fromData())
        node_params['tau'] <- round(1/summary(fit)$sigma, getOption("Hyde_maxDigits"))
      else if (node_params['tau'] == fromFormula())
        stop("parameter 'tau' can not be estimated from a formula.")
    }
    
    model_code <- paste0(node_str, " ~ ",
           network$nodeType[[node_str]], 
           "(", 
           paste(node_params, collapse=", "), 
           ")")
  }
  
  #*** Type 'dpois'
  else if (network$nodeType[[node_str]] == "dpois"){
    if (any(c(fromData()) %in% node_params))
      fit <- do.call(network$nodeFitter[[node_str]],
                     c(list(formula = network$nodeFormula[[node_str]],
                            data = if (is.null(network$nodeData[[node_str]])) network$data else network$nodeData[[node_str]]),
                       network$nodeFitterArgs[[node_str]]))
    
    if (node_params['lambda'] %in% c(fromData(), fromFormula())){
      if (node_params['lambda'] == fromData())
        node_params['lambda'] <- writeJagsFormula(fit, network$nodes)
    
      else if (node_params['lambda'] == fromFormula())
        node_params['lambda'] <- rToJags(network$nodeFormula[[node_str]])
    
      node_params['lambda'] <- as.character(as.formula(node_params[['lambda']]))[-(1:2)]
    }
    
    model_code <- paste0(node_str, " ~ ",
                         network$nodeType[[node_str]],
                         "(",
                         paste(node_params, collapse=", "),
                         ")")
  }
  
  else{
    if (any(node_params %in% c(fromData(), fromFormula())))
      stop(paste0("nodeType '", network$nodeType[[node_str]], 
                  "' does not currently support fromData() or fromFormula()"))
    model_code <- paste0(node_str, " ~ ", 
                         network$nodeType[[node_str]],
                         "(",
                         paste(node_params, collapse=", "),
                         ")")
  }
  
  

  #* Return model  
  return(model_code)
}
