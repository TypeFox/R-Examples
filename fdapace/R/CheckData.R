#' Check data format
#' 
#' Check if there are problems with the form and basic structure of the functional data 'y' and the recorded times 't'.
#' 
#' @param y is a n-by-1 list of vectors
#' @param t is a n-by-1 list of vectors


CheckData = function(y,t){
  
  if(!is.list(y)){
    stop('y should be list \n')
  }
  if(!is.list(t)){
    stop('t should be list \n')
  }
 
 if( length(t) != length(y)){
    stop('t and y should have the same length \n')
  }
 
  #if(any(is.na(unlist(y)))){
  #  stop('y cannot contain NA/NaN entries\n')
  #}
  #if(any(is.na(unlist(t)))){   
  #  stop('t cannot contain NA/NaN entries\n')
  #} 

  ni_y = unlist(lapply(y,function(x) sum(!is.na(x))))
  if(all(ni_y == 1)){  
    stop("FPCA is aborted because the data do not contain repeated measurements in y!\n"); 
  }
  ni_tt = unlist(lapply(t,function(x) sum(!is.na(x))))
  if(all(ni_tt == 1)){  
    stop("FPCA is aborted because the data do not contain repeated measurements in t!\n"); 
  }   
  if( !all(unlist(lapply(y,function(x) class(x) %in% c('integer', 'numeric') ) ) ) ){
        stop("FPCA is aborted because 'y' members are not all of class numeric! Try  \"lapply(y,function(x) class(x))\" to see the current classes. \n");
  }
 if( !all(unlist(lapply(t,function(x) class(x) %in% c('integer', 'numeric'))) ) ){
        stop("FPCA is aborted because 't' members are not all of class numeric! Try  \"lapply(t,function(x) class(x))\" to see the current classes. \n");
  }

 if(any( unlist( lapply(t, function(x) length(x) != length(unique(x))))) ){
        stop("FPCA is aborted because within-subject 't' members have duplicated values.  Try  \"which( unlist( lapply(t, function(x) length(x) != length(unique(x)))))\" to see potentially problematic entries. \n");
  }

}

