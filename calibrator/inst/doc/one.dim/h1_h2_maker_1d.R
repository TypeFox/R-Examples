# This file is intended to be called by calex_1d.R.  It defines basis
# functions for the code, H1.1d(), and the model inadequacy, H2.1d():

h1.1d <- function(x){
out <- c(1,x)
names(out) <- c("const" , "x", "A")
return(out)
}

H1.1d <- function(D1){

    if (is.vector(D1)) {
        D1 <- t(D1)
    }
    out <- t(apply(D1, 1, h1.1d))
    colnames(out) <- names(h1.1d(D1[1,,drop=TRUE]))
    return(out)
}

h2.1d <- function(x){
  ##  out <- c(sin(x[1]*4),x[1]^3)
  f <- function(x){(x-0.5)^2-1/6}
  g <- function(x){sin(x*pi*2)}
  
  out <- c(f(x[1]),g(x[1]))
  names(out) <- c("parabolic","periodic")

  return(out)
}

H2.1d <- 
function (D2) 
{
    if (is.vector(D2)) {
        D2 <- t(D2)
    }
    out <- t(apply(D2, 1, h2.1d))
    colnames(out) <- names(h2.1d(D2[1, , drop = TRUE]))
    return(out)
}

