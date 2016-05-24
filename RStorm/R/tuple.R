
Tuple <- function(x, ...){
  if(!is.data.frame(x)){
    stop("Tuples should be a single row dataframe")
  }
  if(nrow(x) != 1){
    stop("Tuples should be a single row dataframe")
  }
  class(x) <- "Tuple"
  x
}

is.Tuple <- function(x){
  ifelse( class(x) == "Tuple", TRUE, FALSE )
}
