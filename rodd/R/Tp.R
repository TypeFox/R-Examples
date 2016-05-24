Tp <- function(x, w, p, eta, theta.fix, theta.var)
{
    res <- 0
    for(i in 1:length(eta))
        for(j in 1:length(eta))
            if(p[i,j] != 0)
                res <- res + p[i,j] * Tfs(x, w, eta[[i]], eta[[j]], theta.fix[[i]], theta.var[[i,j]])
    res
}
