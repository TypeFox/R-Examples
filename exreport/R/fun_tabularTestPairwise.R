#' Display pairwise information about a multiple test between the methods
#'
#' This function obtain a pairwise table comparing the methods among themselves
#' for the specified metrics. It takes an testMultiplePairwise object as an
#' input.
#'
#' @export
#' @param ph The input testMultiplePairwise object
#' @param value Indicates the metric to be displayed ("pvalue", "wtl")
#' @param charForNAs Indicates the character included when there is not 
#' comparison available
#' @return An extabular object
#' 
#' 
#' @examples
#' # First we create an experiment from the wekaExperiment problem and prepare
#' # it to apply the test:
#' experiment <- expCreate(wekaExperiment, name="test", parameter="fold")
#' experiment <- expReduce(experiment, "fold", mean)
#' experiment <- expInstantiate(experiment, removeUnary=TRUE)
#' 
#' # Then we perform a a testMultiplePairwise test procedure
#' test <- testMultiplePairwise(experiment, "accuracy", "max")
#' 
#' # Different tables can be obtained by using a range of metrics
#' tabularTestPairwise(test, "pvalue")
#' 
#' tabularTestPairwise(test, "wtl")
#' 
tabularTestPairwise <- function(ph, value="pvalue", charForNAs = "-")  {

  if( !is.testMultiplePairwise(ph) )
    stop(.callErrorMessage("wrongParameterError", "ph", "testMultiplePairwise"))
  
  if( !(value %in% c("pvalue","wtl")))
    stop(.callErrorMessage("parameterInRangeError", "value", "from the list [pvalue,wtl]"))
  
  # Compute the values by calling the proper metric function
  metric <- do.call(.mapId2Name[[value]],list(ph))
  
  # Generate data.frame from method names and metric value
  mValues <- metric$data
  mFormat  <- metric$format
  # First if metric has more than one column make a concatenated text version.
  # As we have to interpret the format of each column in order to build only 1 column,
  # we assume that if the all the elements from the row r share the presentation (i.e. bf) 
  if (ncol(mValues) > 1)
  {
    aux = data.frame(values=rep("", nrow(mValues)), stringsAsFactors = FALSE)
    for (i in 1:nrow(mValues)){
      values <- unlist(mValues[i,])
      format <- unlist(mFormat[i,])
      aux[i,] <- paste(sprintf(format,values),collapse = "/")
    }
    
    mValues <- aux
    mFormat <- data.frame(rep("%s",nrow(mValues)))
    colnames(mFormat) <- colnames(mValues)
  }

  d <- cbind(ph$names, mValues)
  dFormat <- cbind(ph$names, mFormat)
  
  # Generate tabular form
  t <- reshape2::dcast(d, method1~method2, value.var=names(d)[3])
  tFormat <- reshape2::dcast(dFormat, method1~method2, value.var=names(d)[3])
  
  # Compute proper order of the matrix:
  means <- rowMeans(ph$friedman$ranks)
  means <- means[order(means)]
  
  # Replace column 1 with names
  method <- t[,1]
  t <- t[,-1]
  rownames(t) <- method
  tFormat <- tFormat[,-1]
  rownames(tFormat) <- method
  
  # Reorder columns and rows
  tab <- t[names(means)[1:length(means)-1],names(means)[2:length(means)]]
  tabFormat <- tFormat[names(means)[1:length(means)-1],names(means)[2:length(means)]]
  
  # Regenerate first row as rownames
  tab <- cbind(method=rownames(tab), tab)
  rownames(tab) <- NULL
  tabFormat <- cbind(method=rep("%s",nrow(tabFormat)), tabFormat)
  rownames(tabFormat) <- NULL
  
  
  # For convenience, we deal with "method" and "method1" and "method2" names for columns
  # of names data.frame. But now, to generate the table, we use the real value (definded by the user).
  methodName <- ph$experiment$method
  colnames(tab)[1] <- methodName
  
  # Customize table if pvalue is present:
  if (value=="pvalue") {
    title <- sprintf("Summary of pairwise post-hoc test for output \"%s\"", 
                     ph$tags$target)
    tableType = "phtest"
  } else {
    title <- sprintf("Method comparison for output \"%s\"", 
                     ph$tags$target)
    tableType = "plain"
  }
  
  res <- .exTabular(tables = list("testPairwise"=tab), 
                    formats=list("testPairwise"=tabFormat), 
                    tableSplit=1, 
                    tableType=tableType, 
                    title = title, 
                    tags = ph$tags)
  res
}
