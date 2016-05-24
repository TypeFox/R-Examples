rdirichlet=function (n, par) 
{
    k = length(par)
    z = array(0, dim = c(n, k))
    s = array(0, dim = c(n, 1))
    for (i in 1:k) {
        z[, i] = rgamma(n, shape = par[i])
        s = s + z[, i]
    }
    for (i in 1:k) {
        z[, i] = z[, i]/s
    }
    return(z)
}
