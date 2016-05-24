mp8.ll <-
function(theta, nrisk, ndeath, age = c(1e-05, 1, seq(5, 85, 5))) 
{
    lx <- nrisk
    dx <- ndeath
    p.hat <- mod8p(theta = theta, x = age)
    ll <- ll.binom(x = dx, n = lx, p = p.hat)
    return(ll)
}

