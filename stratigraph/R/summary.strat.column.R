"summary.strat.column" <-
function(object, ...){
  if(is.strat.column(object)){
    str(object)
  }else{
    stop('this object is not of class strat.column')
  }
} # End of function