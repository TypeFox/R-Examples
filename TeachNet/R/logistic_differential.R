logistic.differential <- function(x) {
  f <- x
  f <- vapply(f, function(m) if(is.infinite(exp(-m))){return(0)}else{return(exp(-m)/(1+exp(-m))^2)},1)
  return(f)
}
