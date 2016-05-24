numberof <- function(x,f) { 
  sum(rep(1,length(x))[f(x)])
}
