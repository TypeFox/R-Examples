"strat.column" <-
function(counts, ...){
  if(is.strat.column(counts)){
    stop('argument to is.strat.column() is already an strat.column object')
  }
  x <- list(...)
  x$counts <- as.matrix(counts)
  if(is.null(x$taxa)){
    x$taxa <- colnames(counts)
  }
  class(x) <- 'strat.column'
  return(x)
} # End of function