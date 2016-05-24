#' force.balance --- repeatedly applies balance until 
#' sub-tolerance is reached
#' INPUT = network model
#' OUTPUT = balanced model
#' M. Lau 1 Oct 2012
#' ---------------------------------------------------

force.balance <- function(x,tol=5,max.itr=10,method='AVG2'){
  n.itr <- 1 # initiate counter
  while(ssCheck(x)==FALSE & n.itr<max.itr){
    x <- balance(x,method=method,tol=tol)
    n.itr <- n.itr + 1
  }
  if (n.itr>=max.itr){
    warning('Maximum iterations reached.')
  }else{
    return(x)
  }
}
