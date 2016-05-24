#' Check if NaN are present in the data and if yes remove them
#' 
#' Check if there are problems cause by missing values with the form and basic structure of the functional data 'y' and the recorded times 't'.
#' 
#' @param y is a n-by-1 list of vectors
#' @param t is a n-by-1 list of vectors

HandleNumericsAndNAN <- function(y,t){
 
  # Check for the presense of NA and remove them (if they exist) from the two lists in a pairwise manner
  if( any(is.na(unlist(t))) ||  any(is.na(unlist(y))) ){
   
    helperF <- function(x) which(!is.na(unlist(x)))
    L <- list(); for(j in 1:length(y)) L[[j]] = list(y[[j]],t[[j]])
    validIndexes = lapply(L, function(x) intersect(helperF(x[1]), helperF(x[2]) ))

    y = lapply(1:length(y), function(i) y[[i]][validIndexes[[i]]])
    t = lapply(1:length(y), function(i) t[[i]][validIndexes[[i]]])

    if( any(unlist(lapply(y, function(x) length(x) == 0))) ){
       stop('Subjects with only NA values are not allowed.\n')
    }
    
    ni_y = unlist(lapply(y,function(x) sum(!is.na(x))))
    if(all(ni_y == 1)){  
      stop("FPCA is aborted because the data do not contain repeated measurements after removing NA values.\n"); 
    }
  }
  

  
  # Force the data to be list of numeric members
  y <- lapply(y, as.numeric) 
  t <- lapply(t, as.numeric)
  t <- lapply(t, signif, 14)
  return( inputData <- list(y=y, t=t));

}
