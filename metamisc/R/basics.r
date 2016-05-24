logit <- function(x) { if(is.numeric(x))  log(x/(1-x)) else stop("x is not numeric!") }

inv.logit <- function(x) {  if(is.numeric(x)) 1/(1+exp(-x)) else stop("x is not numeric!") }


