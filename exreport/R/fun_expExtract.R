#' Extract statistically equivalent methods from a multiple comparison test
#'
#' This functions generates a new experiment incluing the methods that obtained
#' an equivalent performance with statisticall significance in the multiple
#' comparison test i.e. those whose hypotheses were not rejected
#'
#' @export
#' @param ph A testMultipleControl test object
#' @return an experiment object
#' 
#' 
#' @examples
#' # First we create an experiment from the wekaExperiment problem and prepare
#' # it to apply the test:
#' experiment <- expCreate(wekaExperiment, name="test", parameter="fold")
#' experiment <- expReduce(experiment, "fold", mean)
#' experiment <- expInstantiate(experiment, removeUnary=TRUE)
#' 
#' # Then we perform a testMultiplePairwise test procedure
#' test <- testMultipleControl(experiment, "trainingTime", "min")
#' 
#' expExtract(test)
#' 
expExtract <- function(ph){
  # Extracts the experiment from a post hoc control test object. It is a subset of the experiment that was used
  # by the test, remaining the experiments whose p-value is equal or greater than the alpha used.
  #
  # Args:
  #   ph:          the testMultipleControl object
  #
  # Returns:
  #   a experiment object to be used in the toolbox
  
  # PARAMETER VALIDATION:
  # Check if parameters are correct
  if (!is.testMultipleControl(ph))
    stop(.callErrorMessage("wrongParameterError", "ph", "testMultipleControl"))
  
  e <- ph$experiment
  remaining <- ph$names[ph$pvalue>=ph$tags$alpha | is.na(ph$pvalue),]
  
  # As well as in expInstantiate, we mix the parameters and the method (to compare names)
  e$data <- e$data[interaction(e$data[,c(e$method,e$parameters)],sep = ",",drop=TRUE) %in% remaining,,drop=FALSE]
  
  #Append this operation in the historic
  phDescription <- sprintf("%s post-hoc test with %s p-value adjustment for output %s", ph$tags$scope, ph$tags$method, ph$tags$target)
  e$historic <- c(e$historic, list(paste("Subset of experiment based on the results of the", phDescription,sep=" ")))
  
  e
}