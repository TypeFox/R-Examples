n.knots=function(n, cutoff=35, rate=0.2){
    as.integer(trunc(pmin(n,cutoff+pmax(0,n-cutoff)^rate)))
}
