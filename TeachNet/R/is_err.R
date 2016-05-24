is.err <- function(x){ # test for err.fct
  if(!is.character(x)){stop("The err.fct you entered is not of type character!")}
  if(!is.na(x[2])){stop("The err.fct you entered is a vector!")}
  if(!(x %in% c("sse","ce") )){stop("The err.fct you entered is not valid!")}
}