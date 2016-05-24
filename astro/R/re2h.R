re2h = function(n, re = 1){
    # convert re to scalelength
    bn = qgamma(0.5,2*n)
    h = re/(bn^n)
    return(h)
}

