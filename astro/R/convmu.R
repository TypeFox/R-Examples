convmu = function(r, mu, n, re = 1, rmu = re){
    # convert between surface brightnesses
    bn = qgamma(0.5,2*n)
    mu2 = mu + ((2.5*bn)/log(10))*(((r/re)^(1/n))-((rmu/re)^(1/n)))
    return(mu2)
}

