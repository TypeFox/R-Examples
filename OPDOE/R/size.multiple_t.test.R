size.multiple_t.test=function(a, alpha, beta, delta, sd)
{
    k=power.t.test(sig.level=alpha,power=1-beta,delta=delta,sd=sd)$n
    k=ceiling(k)
    N=k*a
    structure(list(size = k, total=N))
}
