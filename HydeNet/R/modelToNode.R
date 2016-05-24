#' @name modelToNode
#' @export modelToNode
#' @importFrom stats as.formula
#' @importFrom stats terms
#' @importFrom utils tail
#' 
#' @title Convert a Model Object to a Network Node
#' @description In cases where model objects may already be fit and established,
#'   they can be used to generate a network without having to refit the 
#'   models or specify the distributions in \code{setNode}.
#'   
#' @param model A model object
#' @param ... Additional arguments to be passed to other functions.  Currently ignored.
#' @param nodes A vector of node names, usually as \code{network$nodes}
#' 
#' @description For models built with \code{xtabs}, although a data frame may be
#'   passed when building the model, it is not stored in the object.  Thus,
#'   the data used to construct the models are not carried with the node.  However,
#'   the JAGS code is built appropriate from the object and this should be of
#'   little concern.

modelToNode <- function(model, nodes, ...) UseMethod("modelToNode")

#' @rdname modelToNode
#' @export

modelToNode.default <- function(model, nodes, ...){
  fitter <- utils::tail(as.character(as.list(model$call)[[1]]), 1)
  if (!fitter %in% c("lm", "glm", "multinom", "xtabs"))
    stop(paste0("The Hyde package only accepts models built by the following functions:\n",
                "  lm, glm (family=\"binomial\" only), multinom, xtabs"))
}

#' @rdname modelToNode
#' @export

modelToNode.cpt <- function(model, nodes, ...)
{
  if (missing(nodes))
    nodes <- names(dimnames(model))
  
  node_name = utils::tail(names(dimnames(model)), 1)
  list(nodes = node_name,
       parents = if (length(dimnames(model)) == 1)
         NULL
       else 
         names(dimnames(model))[-length(names(dimnames(model)))],
       nodeType = "dcat",
       nodeFormula = stats::as.formula(paste0(utils::tail(names(dimnames(model)), 1),
                                       " ~ ",
                                       paste0(names(dimnames(model))[-length(names(dimnames(model)))],
                                              collapse = " + "))),
       nodeFitter = "cpt",
       nodeFitterArgs = list(data = attributes(model)$model,
                             wt = utils::tail(names(attributes(model)$model), 1)),
       nodeParams = list(p = gsub("[[:print:]]+~", "", 
                                  writeJagsFormula(model, 
                                                   nodes))),
       nodeDecision = FALSE,
       nodeUtility = FALSE,
       fromData = TRUE,
       nodeData = attributes(model)$model,
       factorLevels = levels(attributes(model)$model[[node_name]]),
       nodePolicyValues = levels(attributes(model)$model[[node_name]]),
       nodeModel = model)
}

#' @rdname modelToNode
#' @export

modelToNode.glm <- function(model, nodes, ...){
  if (missing(nodes))
    nodes <- nodeFromFunction(names(attributes(terms(model))$dataClasses))
  list(nodes = as.character(stats::terms(model))[2],
       parents = if (length(names(attributes(stats::terms(model))$dataClasses)[-1]) == 0)
                      NULL
                     else 
                       matchVars(names(attributes(stats::terms(model))$dataClasses)[-1],
                                 nodes),
       nodeType = "dbern",
       nodeFormula = model$call$formula,
       nodeFitter = as.character(model$call)[1],
       nodeFitterArgs = as.list(model$call)[-c(1, which(names(as.list(model$call)) %in% c("formula", "data")))],
       nodeParams = list(p = gsub("[[:print:]]+~", "", 
                                  writeJagsFormula(model, 
                                                   nodes))),
       nodeDecision = FALSE,
       nodeUtility = FALSE,
       fromData = TRUE,
       nodeData = if ("data" %in% names(as.list(model$call)[-c(1, which(names(as.list(model$call)) == "formula"))])){
         if (is.null(model$model)) stats::update(model, model=TRUE)$model
         else model$model
       } else NULL,
       factorLevels = if (is.factor(model$model[[1]])) levels(model$model[[1]]) else NULL,
       nodePolicyValues = if (is.factor(model$model[[1]])) levels(model$model[[1]]) 
                          else unique(model$model[[1]]),
       nodeModel = model)
}

#' @rdname modelToNode
#' @export

modelToNode.lm <- function(model, nodes, ...){
  if (missing(nodes))
    nodes <- nodeFromFunction(names(attributes(stats::terms(model))$dataClasses))
  list(nodes = as.character(stats::terms(model))[2],
       parents = if (length(names(attributes(stats::terms(model))$dataClasses)[-1]) == 0)
                     NULL
                    else 
                       matchVars(names(attributes(stats::terms(model))$dataClasses)[-1],
                                 nodes),
       nodeType = "dnorm",
       nodeFormula = model$call$formula,
       nodeFitter = as.character(model$call)[1],
       nodeFitterArgs = as.list(model$call)[-c(1, which(names(as.list(model$call)) %in% c("formula", "data")))],
       nodeParams = list(mu = gsub("[[:print:]]+~", "", 
                                   writeJagsFormula(model,
                                                    nodes)),
                         tau = 1/summary(model)$sigma),
       nodeDecision = FALSE,
       nodeUtility = FALSE,
       fromData = TRUE,
       nodeData = if ("data" %in% names(as.list(model$call)[-c(1, which(names(as.list(model$call)) == "formula"))])){
         if (is.null(model$model)) stats::update(model, model=TRUE)$model
         else model$model
       } else NULL,
       factorLevels = NULL,
       policyValues = stats::quantile(model$model[[1]], probs=c(.25, .5, .75), na.rm=TRUE),
       nodeModel = model)
}

#' @rdname modelToNode
#' @export

modelToNode.multinom <- function(model, nodes, ...){
  if (is.null(model$model)) model <- stats::update(model, model=TRUE)
  if (missing(nodes))
    nodes <- nodeFromFunction(names(attributes(stats::terms(model))$dataClasses))
  list(nodes = as.character(stats::terms(model))[2],
       parents = if (length(names(attributes(stats::terms(model))$dataClasses)[-1]) == 0)
                       NULL
                     else 
                       matchVars(names(attributes(stats::terms(model))$dataClasses)[-1],
                                 nodes),
       nodeType = "dcat",
       nodeFormula = model$call$formula,
       nodeFitter = as.character(model$call)[1],
       nodeFitterArgs = as.list(model$call)[-c(1, which(names(as.list(model$call)) %in% c("formula", "data")))],
       nodeParams = list(pi = gsub("[[:print:]]+~", "", 
                                   writeJagsFormula(model, 
                                                    nodes))),
       nodeDecision = FALSE,
       nodeUtility = FALSE,
       fromData = TRUE,
       nodeData = if (!is.null(model$model)) model$model else NULL,
       factorLevels = if (is.factor(model$model[[1]])) levels(model$model[[1]]) else NULL,
       policyValues = if (is.factor(model$model[[1]])) levels(model$model[[1]])
                      else unique(model$model[[1]]),
       nodeModel = model)
}

#' @rdname modelToNode
#' @export

modelToNode.xtabs <- function(model, nodes, ...){
  if (missing(nodes)) 
    nodes <- nodeFromFunction(names(attributes(model)$dimnames))
  list(nodes = names(attributes(model)$dimnames),
       parents = NULL,
       nodeType = "dcat",
       nodeFormula = attributes(model)$call$formula,
       nodeFitter = as.character(attributes(model)$call)[1],
       nodeFitterArgs = as.list(attributes(model)$call)[-c(1, which(names(as.list(attributes(model)$call)) %in% c("formula", "data")))],
       nodeParams = list(pi = writeJagsFormula(model, nodes)),
       nodeDecision = FALSE,
       nodeUtility = FALSE,
       fromData = FALSE,
       nodeData = NULL,
       factorLevels = names(model),
       policyValues = names(model),
       nodeModel = model)
}

