petrorad = function(n, i = 0.5, re = 1){
    # calculate petrosian radius as multiples of re
    temp = function(n, i, re, res=10, loop=5){
        bn = qgamma(0.5,2*n)
        r = seq(0,10*re,len=res)
        lower = 1
        upper = res
        for(j in 1:loop){
            r = seq(r[lower],r[upper],len=res)
            x = bn*((r/re)^(1/n))
            p = suppressWarnings((2*n*igamma(x,2*n))/((exp(1)^(-x))*(x^(2*n))))
            pmin = which.min(abs((1/p)-i))
            lower = max((pmin-1),1)
            upper = min((pmin+1),length(r))
        }
        return=c(r[pmin])
    }
    prad = as.numeric(Vectorize(temp)(n,i,re=1,res=10,loop=10))
    return(prad)
}

