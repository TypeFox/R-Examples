#' Summarize the result of a multiple comparison statistical test in a table
#'
#' This function builds a table from a testMultiple object, either control or
#' pairwise. The htpotheses are added and compared in the table showing the
#' methods and a range of different metrics than can be added to the table.
#' Also the table shows information about rejected hypotheses.
#'
#' @export
#' @param ph The input testMultiple from which the table is generated
#' @param columns A vector indicating the metrics that will be shown in the table
#' @return an extabular object
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
#' test <- testMultipleControl(experiment, "accuracy", "min")
#' 
#' # Different tables can be obtained by using a range of metrics
#' tabularTestSummary(test, c("pvalue"))
#' 
#' tabularTestSummary(test, c("rank", "pvalue", "wtl"))
#' 
tabularTestSummary <- function(ph, columns=c("pvalue"))  {
  
  if( !is.testMultiple(ph) )
    stop(.callErrorMessage("wrongParameterError", "ph", "testMultiple"))
  
  if( sum(columns %in%  c("pvalue","wtl","rank")) != length(columns) )
    stop(.callErrorMessage("parameterInRangeError", "columns", 
                           "from the list [pvalue,wtl,rank]"))
  
  # df should be a data.frame with the names.
  df <- ph$names
  format=data.frame(matrix("%s", nrow=nrow(df),ncol=ncol(df)))
  colnames(format) <- colnames(df)
  decreasingOrder=c(T)
  # and now we append the metrics
  for ( col in columns )
  {
    m <- do.call(.mapId2Name[[col]],list(ph))
    df <- cbind(df, m$data)
    format <- cbind(format, m$format)
    decreasingOrder <- c(decreasingOrder, m$decreasingOrder)
  }
  # Order in cascade
  correct_order <- .orderCascade(df, decreasing = decreasingOrder)
  df  <- df[correct_order,,drop=FALSE]
  format  <- format[correct_order,,drop=FALSE]
  
  # For convenience, we deal with "method" and "method1" and "method2" names for columns
  # of names data.frame. But now, to generate the table, we use the real value (definded by the user).
  methodName <- ph$experiment$method
  if(ncol(ph$names)>1){
    colnames(df)[1:ncol(ph$names)] <- paste(methodName,1:ncol(ph$names))
  }
  else{
    colnames(df)[1] <- methodName
  }
  
  # Customize table if pvalue is present:
  if ("pvalue"%in%columns) {
    title <- sprintf("Summary of %s post-hoc test for output \"%s\"", ph$tags$scope, ph$tags$target)
    tableType = "phtest"
  } else {
    title <- sprintf("Method comparison for output \"%s\"", ph$tags$target)
    tableType = "plain"
  }
  
  tab <- .exTabular(tables = list("testMultiple"=df), 
                    formats =list("testMultiple"=format), 
                    tableSplit=1, 
                    tableType=tableType, 
                    title = title, 
                    alias = "TestSummaryTable",
                    tags = ph$tags)
  tab
}