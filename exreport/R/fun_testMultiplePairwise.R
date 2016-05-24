#' Multiple Comparison Statistical Test (Friedman + Pairwise Shaffer PostHoc)
#'
#' This function perfoms a multiple comparison statistical test for the given
#' experiment. First of all it performs a Friedman Test over all methods. In the 
#' case this test is rejected, meaning that significant differences are present
#' among the methods a post-hoc test is then executed. For that, each pair of
#' methods are compared between each other, and finally a Shaffer familywise 
#' error correction is applied to the resulting p-values.
#'
#' @export
#' @param e Input experiment
#' @param output The output for which the tet will be performed.
#' @param rankOrder The optimization strategy, can be either maximizing "max"
#' or minimizing "min" the target output variable.
#' @param alpha The significance level used for the whole testing procedure.
#' @return an testMultiplePairwise object
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
#' test <- testMultiplePairwise(experiment, "accuracy", "max")
#' 
#' summary(test)
#' 
testMultiplePairwise <- function(e , output, rankOrder = "max", alpha = 0.05) {
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
    warning( sprintf("Friedman test was not rejected with p-value %.4e >= %.4e alpha, computed p-values will not be significant.", f$pvalue, alpha))
  
  # Apply a posthoc test pairwise general to a testFriedman object.
  # The parameters are the same as in postHocTestPairwise function.
  generalPairwise <- .statisticsGeneralPairwise(f)
  
  # Adjust p-values:
  #NAs are ignored when adjusting p-values
  adjustedPvalues <- .statisticsShafferAdjustemnt(generalPairwise$pvalue)
  
  #The method names
  methodNames <- generalPairwise$names
  
  ph <- .testMultiplePairwise(
    names      = methodNames,
    pvalues    = adjustedPvalues,
    friedman   = f,
    experiment = e,
    alpha      = alpha,
    target     = output,
    scope      = "Pairwise",
    method     = "Shaffer",
    tags       = e$tags
  )
  
  ph
}