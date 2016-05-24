# -------------------------------------------------
.phMetricPvalues <- function(x){
  # Extracts the pvalues from a testMultiple object.
  #
  # Args:
  #   x:          the testMultiple object.
  #
  # Returns:
  #   a numeric array of p-values.
  
  # PARAMETER VALIDATION:
  #Check that variable x is a testMultiple object  
  if( !is.testMultiple(x) )
    stop(.callErrorMessage("wrongParameterError", "x", "testMultiple"))
  
  # There is always a correspondence between x$pvalues an x$names 
  pv     <- data.frame(pvalue = x$pvalues)
  format <- data.frame(pvalue = ifelse(pv$pvalue>=x$tags$alpha,"b%.4e","%.4e"))
  
  metric <- .phMetric(pv, format=format, decreasingOrder = TRUE)
  
  metric
}

.phMetricWTL <- function(x){
  # Extracts the pvalues from a testMultiple object.
  #
  # Args:
  #   x:          the testMultiple object.
  #
  # Returns:
  #   a data.frame with the win-tie-loss metric
  
  # PARAMETER VALIDATION:
  #Check that variable x is a testMultiple object  
  if( !is.testMultiple(x) )
    stop(.callErrorMessage("wrongParameterError", "x", "testMultiple"))
  
  ranks    <- x$friedman$ranks
  names    <- x$names
  
  k <- dim(names)[1]
  
  if( is.testMultipleControl(x) )
    strNames <- data.frame(method1 = rep(as.character(x$control), nrow(names)),
                           method2 = as.character(names[,"method"]),
                           stringsAsFactors = FALSE)
  else
    strNames <- data.frame(method1 = as.character(names[,"method1"]),
                           method2 = as.character(names[,"method2"]),
                           stringsAsFactors = FALSE)
  
  win     <- c()
  tie     <- c()
  loss    <- c()
  for (i in 1:k){
    method1 <- strNames[i,"method1"]
    method2 <- strNames[i,"method2"]
    if(method1==method2){
      win[i]        <- NA
      tie[i]        <- NA
      loss[i]       <- NA
    }
    else{
      win[i]        <- sum(ranks[method2,] > ranks[method1,])
      tie[i]        <- sum(ranks[method2,] == ranks[method1,])
      loss[i]       <- sum(ranks[method2,] < ranks[method1,])
    }
  }
  
  dat <- data.frame(win = win, tie = tie, loss = loss)
  format <- data.frame(matrix("%d",nrow = k, ncol=3),stringsAsFactors = FALSE)
  metric <- .phMetric(dat, format=format, decreasingOrder = c(TRUE,TRUE,FALSE))
  metric
}

.phMetricRanks <- function(x){
  # Extracts the ranks from a testMultiple object.
  #
  # Args:
  #   x:          the testMultiple object.
  #
  # Returns:
  #   a numeric array with the ranking
  
  # PARAMETER VALIDATION:
  #Check that variable x is a testMultiple
  if( !is.testMultipleControl(x))
    stop(.callErrorMessage("wrongParameterError", "x", "testMultipleControl"))
  
  # There is always a correspondence between x$ranks an x$names 
  dat <- data.frame(rank = rowMeans(x$friedman$ranks))
  format <- data.frame(matrix("%.2f",nrow = nrow(dat), ncol=1),stringsAsFactors = FALSE)
  res <- .phMetric(dat, format=format, decreasingOrder = FALSE)
  res
}