indepmetrop=function (logpost, proposal, start, m, ...) 
{
    logmultinorm = function(x, m, v) {
        return(-0.5 * t(x - m) %*% solve(v) %*% (x - m))
    }
    pb = length(start)
    Mpar = array(0, c(m, pb))
    mu = matrix(proposal$mu)
    if(diff(dim(mu))>0) mu=t(mu)
    v = proposal$var
    a = chol(v)
    f0 = logpost(start, ...)
    th0 = matrix(t(start))
    accept = 0
    for (i in 1:m) {
        th1 = mu + t(a) %*% array(rnorm(pb), c(pb, 1))
        f1 = logpost(t(th1), ...)
        R = exp(logmultinorm(th0, mu, v) - logmultinorm(th1, 
            mu, v) + f1 - f0)
        u = runif(1) < R
        if (u == 1) {
            th0 = th1
            f0 = f1
        }
        Mpar[i, ] = th0
        accept = accept + u
    }
    accept = accept/m
    return(list(par = Mpar, accept = accept))
}
