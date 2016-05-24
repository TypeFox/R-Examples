#' @name HydeNetwork
#' @export HydeNetwork
#' @importFrom graph nodes
#' @importFrom gRbase dag
#' @importFrom gRbase graphNEL2adjMAT
#' @importFrom stats as.formula
#' 
#' @title Define a Probablistic Graphical Network
#' @description Using either a directed acyclic graph (DAG) or a list of models, 
#'   define a probabilistic
#'   graphical network to serve as the basis of building a model.  
#'   
#' @param nodes Either a formula that defines the network as passed to 
#'   \code{gRbase::dag} or a list of model objects.
#' @param data A data frame with the data for estimating node parameters.
#' @param ... additional arguments to \code{gRbase::dag}.
#' 
#' @details The DAG becomes only one element of the object returned by 
#'   \code{HydeNetwork}. The dag object is used to extract the node names
#'   and a list of parents for each node.  These will be used to help quantify
#'   the relationships.
#'   
#'   When given a formula, the relationships are defined, but are not quantified
#'   until \code{writeNetworkModel} is called.
#'   
#'   When a list of models is given, rather than refitting models when 
#'   \code{writeNetworkModel} is called, the quantified relationships are 
#'   placed into the object.
#'   
#' @return
#' Returns an object of class \code{HydeNetwork}. The object is really just a 
#' list with the following components:
#' \itemize{
#'   \item \code{nodes} a vector of node names
#'   \item \code{parents} a named list with each element being a vector of parents
#'       for the node named.
#'   \item \code{nodeType} a named list with each element specifying the JAGS 
#'       distribution type.
#'   \item \code{nodeFormula} a named list with the formulae specifying the 
#'       relationships between nodes.
#'   \item \code{nodeFitter} a named list giving the fitting function for each
#'       node.
#'   \item \code{nodeFitterArgs} A named list with additional arguments to be passed
#'       to fitter functions.
#'   \item \code{nodeParams} A named list.  Each element is a vector of parameters 
#'       that will be expected by JAGS.
#'   \item \code{fromData} A named list with the logical value of wheather parameters
#'       should be estimated from the data.
#'   \item \code{nodeData} A named list with the data for each node.  If a node's 
#'       entry in \code{fromData} is \code{TRUE} and \code{nodeData} is \code{NULL},
#'       it will look to the \code{data} attribute instead.
#'   \item \code{factorLevels} If the vector associated with the node is a factor 
#'       (or character), the levels of the factor are stored here.  Although it
#'       may seem redundant, it allows factor levels to be specified in cases
#'       where the node is not define with data.  If data are provided to the 
#'       node, this element is determined from the data and cannot be 
#'       manually overwritten.
#'   \item \code{nodeModel} A list of model objects.  This is a storing place for 
#'       models that have already been fit so that they don't have to be refit 
#'       again.
#'   \item \code{nodeDecision} A named list of logical flags for whether the node is
#'       a decision node or not.
#'   \item \code{nodeUtility} A named list of logical flags for whether the node is
#'       a utility node or not.
#'   \item \code{dag} The dag object returned by \code{gRbase}.  Most of the plotting
#'       utilities will be based on this element.
#'   \item \code{data} A common data frame for nodes that do not have their own unique
#'       data source.
#'   \item \code{network_formula} The original formula passed to \code{gRbase::dag}
#'       to construct the model.
#'  }
#'  
#'  @note These objects can get pretty large.  In versions of R earlier than 3.2, 
#'    it can take a while to print the large network objects if you simply type
#'    the object name into the console.  It is recommended that you always 
#'    explicitly invoke the `print` function (ie, \code{print(Net)} instead
#'    of just \code{Net}) to save yourself some valuable time.
#'   
#' @author Jarrod Dalton and Benjamin Nutter
#' @examples
#' #* Formula Input
#' Net <- HydeNetwork(~ wells + 
#'                      pe | wells + 
#'                      d.dimer | pregnant*pe + 
#'                      angio | pe + 
#'                      treat | d.dimer*angio + 
#'                      death | pe*treat,
#'                      data = PE) 
#' print(Net)
#' 
#' #* Model Input
#' g1 <- lm(wells ~ 1, data=PE)
#' g2 <- glm(pe ~ wells, data=PE, family="binomial")
#' g3 <- lm(d.dimer ~ pe + pregnant, data=PE)
#' g4 <- xtabs(~ pregnant, data=PE)
#' g5 <- glm(angio ~ pe, data=PE, family="binomial")
#' g6 <- glm(treat ~ d.dimer + angio, data=PE, family="binomial")
#' g7 <- glm(death ~ pe + treat, data=PE, family="binomial")
#' 
#' bagOfModels <- list(g1,g2,g3,g4,g5,g6,g7)
#' 
#' bagNet <- HydeNetwork(bagOfModels)
#' print(bagNet)
#' 


HydeNetwork <- function(nodes, ...) UseMethod("HydeNetwork")

#' @rdname HydeNetwork
#' @export


HydeNetwork.formula <- function(nodes, data=NULL, ...){
  #* Build the DAG object
  network <- gRbase::dag(nodes) 
  
  #* Node names
  node_names <- graph::nodes(network)
  
  #* Parents
  parents <- HydeNetwork_parents(network)
  names(parents) <- node_names
  
  #* fromData
  #* returns TRUE if the node and its parents are in 'data'
  #* returns FALSE if any node or parent is missing from 'data'
  fromData <- lapply(node_names, HydeNetwork_fromData, data, parents)
  names(fromData) <- node_names
  
  #* nodeFormula
  nodeFormula <- lapply(seq_along(parents), HydeNetwork_nodeFormula, parents, data, fromData)
  names(nodeFormula) <- node_names
  
  #* nodeFitter
  #* returns 'lm' for continuous variables
  #* returns 'glm' for categorical variables
  #* returns NULL for variables not in the data
  nodeFitter <- lapply(node_names, HydeNetwork_nodeFitter, data, parents)
  names(nodeFitter) <- node_names

  #* nodeTypes
  #* returns 'dcat' if categorical and has no parents
  #* returns 'dnorm' otherwise
  nodeType <- lapply(node_names, HydeNetwork_nodeType, data, parents, nodeFitter)
  names(nodeType) <- node_names

  #* nodeParameters
  data(jagsDists, envir=environment())
  nodeParams <- lapply(node_names, HydeNetwork_nodeParams, jagsDists, nodeType, fromData)
  names(nodeParams) <- node_names
  
  #* fitterArgs
  nodeFitterArgs <- lapply(seq_along(node_names), function(x) return(NULL))
  if (any(sapply(nodeFitter, nodeFitter_is_glm)))
    nodeFitterArgs[which(sapply(nodeFitter, function(x) x == "glm"))] <- list(family='binomial')
  names(nodeFitterArgs) <- node_names
  
  nodeData <- lapply(seq_along(node_names), function(x) return(NULL))
  names(nodeData) <- node_names
  
  nodeModel <- lapply(seq_along(node_names), function(x) return(NULL))
  names(nodeModel) <- node_names  
  
  nodeDecision <- lapply(seq_along(node_names), function(x) return(FALSE))
  names(nodeDecision) <- node_names 
  
  nodeUtility <- lapply(seq_along(node_names), function(x) return(FALSE))
  names(nodeUtility) <- node_names
  
  factorLevels <- lapply(seq_along(node_names), function(x) return(NULL))
  names(factorLevels) <- node_names
  if (!is.null(data)){
    factor_vars <- names(data)[vapply(data, is.factor, logical(1))]
    factorLevels[factor_vars] <- 
      lapply(data[, factor_vars, drop = FALSE],
             levels)
  }
  
  policyValues <- lapply(node_names, HydeNet_defaultPolicyValues, data)
  names(policyValues) <- node_names

  #* Define the HydeNetwork object
  network <- list(nodes = node_names, parents=parents, nodeType=nodeType,
                  nodeFormula=nodeFormula,
                  nodeFitter=nodeFitter, nodeFitterArgs=nodeFitterArgs,
                  nodeParams=nodeParams, 
                  fromData=fromData, 
                  nodeData = nodeData,
                  factorLevels = factorLevels,
                  nodeModel = nodeModel,
                  nodeDecision = nodeDecision,
                  nodePolicyValues = policyValues,
                  nodeUtility = nodeUtility,
                  dag=network)
  
  
  network$data <- if (!is.null(data)) data else NULL
  network$network_formula <- nodes
  class(network) <- c("HydeNetwork")
  return(network)
}

#' @rdname HydeNetwork
#' @export

HydeNetwork.list <- function(nodes, ...){
  #* convert models to nodes
  Attrs <- lapply(nodes, modelToNode)
  
  #* assign names to list elements
  for(i in 1:length(Attrs)){
    names(Attrs)[i] <- Attrs[[i]]$nodes  
  }
  
  #* Generate the DAG formula and build the network
  dag.form <- sapply(Attrs, 
                     function(x) paste0(x$nodes, 
                                        if (!is.null(x$parents)) " | " else "", 
                                        paste(x$parents, collapse=" * ")))
  dag.form <- paste0("~ ", paste(dag.form, collapse = " + "))
  network <- HydeNetwork(stats::as.formula(dag.form))
  
  #* Reassign parameters from the models
  for (i in names(Attrs)){
    network$nodeType[[i]] <- Attrs[[i]]$nodeType
    network$nodeFormula[[i]] <- Attrs[[i]]$nodeFormula
    network$nodeFitter[[i]] <- Attrs[[i]]$nodeFitter
    network$nodeFitterArgs[[i]] <- Attrs[[i]]$nodeFitterArgs
    network$nodeParams[[i]] <- Attrs[[i]]$nodeParams
    network$nodeData[[i]] <- Attrs[[i]]$nodeData
    network$nodeModel[[i]] <- Attrs[[i]]$nodeModel
    network$nodeDecision[[i]] <- Attrs[[i]]$nodeDecision
    network$nodeUtility[[i]] <- Attrs[[i]]$nodeUtility
    network$fromData[[i]] <- TRUE
    network$factorLevels[[i]] <- Attrs[[i]]$factorLevels
    network$nodePolicyValues[[i]] <- Attrs[[i]]$policyValues
  }
  
  return(network)
}

#** Utility Functions ***************************
HydeNetwork_parents <- function(network){
  adjMat <- gRbase::graphNEL2adjMAT(network)  
  parents <- lapply(1:ncol(adjMat), function(x) rownames(adjMat)[adjMat[, x] == 1])
  parents <- lapply(parents, function(x) if (length(x) == 0) NULL else x)
  parents
}

HydeNetwork_fromData <- function(node_names, data, parents){
    if (is.null(data)) return(FALSE)
    if (all(c(node_names, parents[[node_names]]) %in% names(data)))
      return(TRUE)
    else return(FALSE)
}

HydeNetwork_nodeFormula <- function(x, parents, data, fromData){
  if (is.null(parents[[x]])){
    if (fromData[[names(parents)[x]]] & !is.numeric(data[, names(parents)[x]]))
      f <- paste("~ ", names(parents)[x])
    else f <- paste(names(parents)[x], "~ 1")
  }
  else f <- paste(names(parents)[x], "~", paste(parents[[x]], collapse=" + "))
  return(stats::as.formula(f))
}

HydeNetwork_nodeFitter <- function(node_name, data, parents){
  if (is.null(data)) return(NULL)
  if (!node_name %in% names(data)) return(NULL)
  else if (is.numeric(data[, node_name])) return("lm")
  else if (is.factor(data[, node_name]) & is.null(parents[[node_name]])) return("xtabs")
  else if (is.factor(data[, node_name]) & 
           all(vapply(parents[[node_name]], function(p) is.factor(data[, p]), logical(1))))
    return("cpt")
  else if (is.factor(data[, node_name]) & nlevels(data[, node_name]) == 2) return("glm")
  else if (is.factor(data[, node_name]) & nlevels(data[, node_name]) > 2) return("multinom")
  else return("glm")
}

HydeNetwork_nodeType <- function(node_name, data, parents, nodeFitter){
  if (is.null(data)) return('dnorm')
  if (node_name %in% names(data)){
    if ((is.null(parents[[node_name]]) && 
         !is.numeric(data[, node_name])) || 
        (!is.null(parents[[node_name]]) && 
         !is.numeric(data[, node_name]) && 
         nlevels(data[, node_name]) > 2))
      return('dcat')
    else if (nodeFitter[[node_name]] == "cpt") return('dcat')
    else if ((is.null(parents[[node_name]]) && 
              !is.numeric(data[, node_name])) || 
             (!is.null(parents[[node_name]]) && 
              !is.numeric(data[, node_name]) && 
              nlevels(data[, node_name]) == 2))
      return('dbern')
    else return('dnorm')
  }
  else return('dnorm')
}  

HydeNetwork_nodeParams <- function(node_name, jagsDists, nodeType, fromData){
  parm <- jagsDists$Parameters[jagsDists$FnName == nodeType[[node_name]]]
  if (fromData[[node_name]]) 
    parm <- paste0("c(",
                   paste(parm, "fromData()", sep="=", collapse=", "),
                   ")")
  else 
    parm <- paste0("c(",
                   paste(parm, "'Unspecified'", sep="=", collapse=", "),
                   ")")
  return(eval(parse(text=parm)))
}

HydeNet_defaultPolicyValues <- function(node_name, data){
  if (is.null(data)) return(NULL)
  if (!node_name %in% names(data)) return(NULL)
  else {
    if (is.numeric(data[[node_name]]))
      return(stats::quantile(data[[node_name]], probs = c(.25, .50, .75), na.rm=TRUE))
    else if (is.factor(data[[node_name]]))
      return(levels(data[[node_name]]))
    else
      return(unique(data[[node_name]]))
  }
}

nodeFitter_is_glm <- function(fitter){
  if (is.null(fitter)) FALSE else fitter == "glm"
}
