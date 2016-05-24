CumInt <-
function(x,int,...){
  fn <- function(x)
    integrate(function(u)sapply(u,int),subdivisions=1000,lower=0,upper=x,...)$val
  sapply(x,fn)
}

