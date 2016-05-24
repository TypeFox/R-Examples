"is.strat.column" <-
function(x, ...){
  if(class(x) == 'strat.column'){
    if(length(dim(x$counts)) == 2 && is.numeric(x$counts)){
      return(TRUE)
    }else{
      cat('argument to is.strat.column() has class attribute \'strat.column\', but x$counts is not a 2-dimensional numeric array\n')
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
} # End of function