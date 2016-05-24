"as.strat.column" <-
function(x){
  class(x) <- 'strat.column'
  if(is.strat.column(x)){
    return(x)
  }else{
    stop('x seems to lack a valid count matrix and cannot be coerced to a strat.column')
  }
} # End of function