is.learn <- function(x) { # test for learning.rate
  if(!is.numeric(x)){stop("The learning.rate you entered is not numeric!")}
  if(!is.na(x[2])){stop("The learning.rate you entered is a vector!")}
  if(!(x>0)){stop("The learning.rate you entered is zero or less!")}
  if(!(x<1)){stop("The learning.rate you entered is one or more!")}
}