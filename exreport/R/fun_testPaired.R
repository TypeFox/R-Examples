#' Paired Wilcoxon statistical test
#'
#' This function performs a Wilcoxon paired test to compare the methods of an
#' experiment consisting exactly on two of them. If more methods are present,
#' then a multiple comparison test must be applied.
#'
#' @export
#' @param e Input experiment
#' @param output The output for which the tet will be performed.
#' @param rankOrder The optimization strategy, can be either maximizing "max"
#' or minimizing "min" the target output variable.
#' @param alpha The significance level used for the whole testing procedure.
#' @return a testPaired object
#' 
#' 
#' @examples
#' # First we create an experiment from the wekaExperiment problem and prepare
#' # it to apply the test, we must subset it to only two methods:
#' experiment <- expCreate(wekaExperiment, name="test", parameter="fold")
#' experiment <- expSubset(experiment, list(method = c("J48", "NaiveBayes")))
#' experiment <- expSubset(experiment, list(featureSelection = c("no")))
#' experiment <- expReduce(experiment, "fold", mean)
#' experiment <- expInstantiate(experiment, removeUnary=TRUE)
#' 
#' # Then we perform a Wilcoxon test procedure
#' test <- testPaired(experiment, "accuracy", "max")
#' 
#' summary(test)
#' 
testPaired <- function(e , output, rankOrder = "max", alpha = 0.05){
  # PARAMETER VALIDATION:
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
  
  #Extracts the problem with only the interested output variable
  data <- e$data[,c(e$method,e$problem,output)]
  #Now transform the problem into a method-versus-problem one.
  data <- reshape2::dcast(data, paste(e$method,"~",e$problem),value.var = output)
  #First column is now the rownames (and then is removed).
  rownames(data) <- data[,1]
  #The parameter drop = FALSE avoid to auto-convert the data.frame into
  #a single array when only one column data.frame is the result of the subseting.
  data <- data[,-1, drop=FALSE]
  
  # Check if it only exists two methods
  if( nrow(data) !=2 )
    stop(.callErrorMessage("wilcoxonNMethodsError"))
  
  
  if (rankOrder == "min" ){
    rankModifier <- 1
  }
  else{
    rankModifier <- -1
  }
  
  ranks                 <- sapply(data*rankModifier, FUN=rank)
  rownames(ranks)       <- rownames(data)
  friedmanRanks         <- rowMeans(ranks)
  
  if(friedmanRanks[1]<friedmanRanks[2]){
    bestMethod <- as.numeric(data[1,])
    bestMethodName <- rownames(data)[1]
    worstMethod <- as.numeric(data[2,])
    worstMethodName <- rownames(data)[2]
  }
  else{
    bestMethod <- as.numeric(data[2,])
    bestMethodName <- rownames(data)[2]
    worstMethod <- as.numeric(data[1,])
    worstMethodName <- rownames(data)[1]
  }
  
  wilcoxon <- .statisticsWilcoxon(bestMethod, worstMethod)
  
  
  #The method names
  methodNames <- rownames(data)
  
  w <- .testPaired(bestMethod  = bestMethodName, 
                   worstMethod = worstMethodName, 
                   statistic   = wilcoxon$statistic, 
                   pvalue      = wilcoxon$pvalue, 
                   alpha       = alpha,
                   target      = output,
                   objetive    = .mapId2Name[[rankOrder]], 
                   tags        = e$tags)
  
  w
}