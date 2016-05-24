concen = function(a, n){
    # central concentration
    bn = qgamma(0.5,2*n)
    c = (igamma(bn*(a^(1/n)),2*n))/(igamma(bn,2*n))
    return(c)
}

