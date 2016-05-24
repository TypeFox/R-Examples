scmean = function(x, mult = 3, loop = 10){
    m = mean(x, na.rm=TRUE)
    s = sd(x, na.rm=TRUE)
    for(i in 1:(loop-1)){
        x = x[x>(m-(mult*s))]
        x = x[x<(m+(mult*s))]
        m = mean(x, na.rm=TRUE)
        s = sd(x, na.rm=TRUE)
    }
    return(list(m=m,s=s))
}

