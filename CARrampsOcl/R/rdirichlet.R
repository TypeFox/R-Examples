rdirichlet <-
function (n, parms) 
{
# generate n random vectors from dirichlet
# rejection envelope

    l <- length(parms)
    x <- matrix(rgamma(l * n, parms), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
}

