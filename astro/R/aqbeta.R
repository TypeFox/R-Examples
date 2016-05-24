aqbeta = function(k, n, s = c(-1,1), p = NA, corr = TRUE, ...){
    if(is.na(p[1])){
        conf = pnorm(s)
        lab = s
    }else{
        conf = p
        lab = p
    }
    qlo = which(conf<0.5)
    qhi = which(conf>0.5)
    aqb = function(k, n, ...){
        qb = qbeta((conf), k+1, n-k+1, ...)
        if(corr){
            o = k/n
            if(is.nan(o)){
                qb[qlo] = NA
                qb[qhi] = NA
            }else{
                if(any(qb[qlo]>o)){
                    qb[qlo] = o
                }
                if(any(qb[qhi]<o)){
                    qb[qhi] = o
                }
            }
        }
        return(qb)
    }
    err = t(Vectorize(aqb)(k=k,n=n))
    colnames(err) = lab
    return(err)
}

