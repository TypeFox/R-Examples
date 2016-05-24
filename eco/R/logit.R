logit <- function(x)
  return(log(x)-log(1-x))

invlogit <- function(x)
  return(exp(x)/(1+exp(x)))
