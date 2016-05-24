#' Multiple Comparison Statistical Test (Friedman + Control Holm PostHoc)
#'
#' This function perfoms a multiple comparison statistical test for the given
#' experiment. First of all it performs a Friedman Test over all methods. In the 
#' case this test is rejected, meaning that significant differences are present
#' among the methods a post-hoc test is then executed. For that, a comparison
#' using the best method as a control is performed for each other method,
#' finally a Holm familywise error correction is applied to the resulting 
#' p-values.
#'
#' @export
#' @param e Input experiment
#' @param output The output for which the tet will be performed.
#' @param rankOrder The optimization strategy, can be either maximizing "max"
#' or minimizing "min" the target output variable.
#' @param alpha The significance level used for the whole testing procedure.
#' @return an testMultipleControl object
#' 
#' 
#' @examples
#' # First we create an experiment from the wekaExperiment problem and prepare
#' # it to apply the test:
#' experiment <- expCreate(wekaExperiment, name="test", parameter="fold")
#' experiment <- expReduce(experiment, "fold", mean)
#' experiment <- expSubset(experiment, list(featureSelection = "yes"))
#' experiment <- expInstantiate(experiment, removeUnary=TRUE)
#' 
#' # Then we perform a testMultiplePairwise test procedure
#' test <- testMultipleControl(experiment, "accuracy", "max")
#' 
#' summary(test)
#' 
testMultipleControl <- function(e , output, rankOrder = "max", alpha = 0.05) {
  #
  # Parameter validation:
  #
  #Check that variable e is an experiment object  
  if( !is.experiment(e) )
    stop(.callErrorMessage("wrongParameterError", "e", "experiment"))
  
  #Check that output variable exists in the experiment e
  if( !(output %in% colnames(e$data)) ) 
    stop(.callErrorMessage("variableNotPresentError", output))
  
  #Check that rankOrder is one of the two possible values
  if(rankOrder!="max" && rankOrder!="min")
    stop(.callErrorMessage("parameterInRangeError", "rankOrder", "list: [max, min]"))
  
  #Check that variable alpha is a valid significance level for the test
  if( !is.numeric(alpha) || alpha<=0 || alpha>=1  )
    stop(.callErrorMessage("parameterInRangeError", "alpha", "range [0, 1]"))
  
  # Check if instantiation of parameters is needed
  if( length(e$parameters) != 0 )
    stop(.callErrorMessage("requireInstantiationError"))
  
  ##
  ##
  
  # Perform Friedman test:
  f <- .doFriedmanTest(e, output, rankOrder, alpha)
  
  # If the friedman test determines that there is no difference 
  # between methods, we stop here.
  if(f$tags$pvalue >= alpha)
    warning(.callWarningMessage("friedmanTestNotRejected", f$tags$pvalue, alpha))
  
  # Apply a posthoc test control to the testFriedman object.
  controlGeneral <- .statisticsGeneralControl(f)
  
  # Adjust p-values:
  #NAs are ignored when adjusting p-values
  adjustedPvalues <- .statisticsHolmsAdjustemnt(controlGeneral$pvalue)
  
  #The method names
  methodNames <- rownames(f$ranks)
  
  ph <- .testMultipleControl(
                     names      = data.frame(method = methodNames),
                     control    = methodNames[controlGeneral$controlIndex],
                     pvalues    = adjustedPvalues,
                     friedman   = f,
                     experiment = e,
                     alpha      = alpha,
                     target     = output,
                     scope      = "Control",
                     method     = "Holm",
                     tags       = e$tags
                     )
  ph
}