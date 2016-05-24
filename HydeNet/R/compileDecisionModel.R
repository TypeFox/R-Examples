#' @name compileDecisionModel
#' @export
#' 
#' @title Compile JAGS Models to Evaluate the Effect of Decisions in a Network
#' @description Nodes at which a decision can be made, such as the decision to 
#'   test or not test; treat or not treat; or use open or robotic surgery may 
#'   impact the outcome for a subject.  These types of decisions may not be 
#'   truly random and understanding how these decisions may impact downstream
#'   outcomes may be beneficial to making the decision.  Compiling the decision
#'   network permits the network to be evaluated under the conditions of each 
#'   set of decisions separately.
#'   
#' @param network A HydeNet object with decision nodes defined.
#' @param policyMatrix A data frame of policies to apply to decision nodes
#'   for comparing networks under different conditions.  See 
#'   \code{\link{policyMatrix}}.
#' @param ... Additional arguments to pass to \code{jags.model}, excepting
#'   the \code{data} argument.  The \code{data} argument is created by 
#'   \code{compileDecisionModel}, and cannot be passed manually.
#' @param data An optional list of data values to be observed in the nodes.  
#'   It is passed to the \code{data} argument of \code{rjags::jags}. Any
#'   values given in \code{data} will override values provided in 
#'   \code{policyMatrix} with a warning.
#'   
#' @details \code{compileDecisionModel} only accepts nodes of type \code{"dbern"}
#'   (Bernoulli random variable taking either 0 or 1) or \code{"dcat"} 
#'   (categorical variables) as decision nodes.  
#'   When the node is type \code{"dcat"}, the 
#'   decision options are extracted from the JAGS statement returned by 
#'   \code{writeJagsModel}.
#'   
#'   The options for each decision nodes (if there are multiple nodes) are 
#'   combined via \code{expand.grid} to make a table of all possible decisions.
#'   Each row of this table is passed as a list to the \code{data} argument 
#'   of \code{jags.model} (via \code{compileJagsModel}) and a list of JAGS
#'   model objects is returned.  \code{coda.samples} may be run on each of these
#'   models.
#'   
#' @return Returns a list of \code{compiledHydeNetwork} objects.
#' 
#' @author Jarrod Dalton and Benjamin Nutter
#' 
#' @seealso \code{\link{policyMatrix}} \code{\link{compileJagsModel}}
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
#' Net <- setDecisionNodes(Net, treat)
#' plot(Net)
#' 
#' decision1 <- compileDecisionModel(Net)
#'
#' #* An effectively equivalent call as the previous
#' decision2 <- compileDecisionModel(Net, policyMatrix(Net))
#' 
#' #* Using a customized policy matrix
#' #* Note: this is a bit of nonsense--you can't decide if a test is negative
#' #*       or positive, but we'll do this for illustration.
#' custom_policy <- policyMatrix(Net, 
#'                               treat="No", 
#'                               angio = c("Negative", "Positive"))
#' decision3 <- compileDecisionModel(Net, custom_policy) 
#' 
compileDecisionModel <- function(network, policyMatrix = NULL, ..., data = NULL){
  Check <- ArgumentCheck::newArgCheck()
  
  dots <- list(...)
  
  options <- makePolicyMatrix(network, policyMatrix, data, Check)

  ArgumentCheck::finishArgCheck(Check)
  
  cpt_arrays <- makeCptArrays(network)
  
  jags.code <- compileJagsModel(network, ...)

  lapply(options,
         runJagsDecisionModel,
         jags.code,
         cpt_arrays, 
         ...)
  
  
}


#*********** UTILITY FUNCTIONS

makePolicyMatrix <- function(network, policyMatrix, data, argcheck){
  if (is.null(policyMatrix))
  {
    decisionNodes <- names(network$nodeDecision)[sapply(network$nodeDecision, any)]
    
    if (length(decisionNodes) == 0)
      ArgumentCheck::addError(
        msg = "No decision nodes indicated in the network",
        argcheck = argcheck)
    
    if (length(decisionNodes) == 0) break; # The next argument check isn't meaningful
    # when this condition is true.
    
    validDecision <- sapply(network$nodeType[decisionNodes], 
                            function(x) x %in% c("dbern", "dcat", "dbin"))
    
    if (!all(validDecision))
      ArgumentCheck::addError(
        msg = paste0("Only nodes of type 'dcat', and 'dbin' may be decision nodes.\n  ",
                     paste0(names(validDecision)[!validDecision], collapse=", "),
                     " cannot be used as decision nodes."),
        argcheck = argcheck)
    
    if (!all(validDecision)) break; # Avoids defining 'options' when there are invalid decision nodes
    
    options <- lapply(decisionNodes, decisionOptions, network)
    names(options) <- decisionNodes
    
    options <- expand.grid(options, stringsAsFactors=FALSE) 
  }
  else
  {
    if (!is.data.frame(policyMatrix))
      ArgumentCheck::addError(
        msg = "'policyMatrix' must be a data frame",
        argcheck = argcheck)
    if (!is.data.frame(policyMatrix)) break; # avoids defining 'options' when
    # the condition is not satisfied
    options <- policyMatrix
  }
  
  #* This is the part that pushes values from `data` into the 
  #* policy matrix.
  if (!is.null(data)){
    conflicts <- names(data)[names(data) %in% names(options)]
    if (length(conflicts) > 0){
      ArgumentCheck::addWarning(
        msg = paste0("The following variables in 'data' are overriding ",
                     "values in 'policyMatrix': ",
                     paste0(conflicts, collapse = ", ")),
        argcheck = argcheck)
    }
    
    for (i in names(data)){
      options[[i]] <- data[[i]]
    }
    #* Remove duplicated rows
    options <- options[!duplicated(options), , drop = FALSE]
  }
  
  ArgumentCheck::finishArgCheck(argcheck)
  
  options <- lapply(1:nrow(options), 
                    function(i){ 
                      l <- as.list(options[i, , drop=FALSE])
                      nms <- names(l)
                      attributes(l) <- NULL
                      names(l) <- nms
                      l
                    })
  
  return(options)
}

#**



makeCptArrays <- function(network){
  cpt_arrays <- unlist(network$nodeFitter) == "cpt"
  if(any(cpt_arrays)){
    cpt_arrays <- names(cpt_arrays)[cpt_arrays]
    cpt_arrays <- network$nodeModel[cpt_arrays]
    nms <- names(cpt_arrays)
    cpt_arrays <- 
      lapply(names(cpt_arrays),
             function(ca){
               if ("cpt" %in% class(cpt_arrays[[ca]])) return(cpt_arrays[[ca]])
               else{
                 args <- 
                   list(formula = network$nodeFormula[[ca]],
                        data = if (!is.null(network$nodeData[[ca]])) network$nodeData[[ca]]
                        else network$data)
                 if (!is.null(network$nodeFitterArgs[[ca]]))
                   args <- c(args, network$nodeFitterArgs[[ca]])
                 return(do.call("cpt", args))
               }   
             })
    names(cpt_arrays) <- paste0("cpt.", nms)
  } else cpt_arrays = list()
  return(cpt_arrays) 
}

#****

runJagsDecisionModel <- function(o, j, cpt_arrays, ...){
  con <- textConnection(paste0(j$jags$model(),
                                 collapse="\n"))
  o_character <- vapply(o, is.character, logical(1))
  for (i in seq_along(o)){
    if (o_character[i])
      o[[i]] <- j$factorRef[[names(o[i])]]$value[j$factorRef[[names(o[i])]]$label == o[[i]]]
  }
  
  cHN <- list(jags = rjags::jags.model(con,
                                       data = c(o, cpt_arrays),
                                       ...),
              observed = o,
              dag = j$dag,
              factorRef = j$factorRef)
  class(cHN) <- c("compiledHydeNetwork")
  close(con)
  return(cHN)
}  
