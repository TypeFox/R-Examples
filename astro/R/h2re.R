h2re = function(n, h = 1){
    # convert scalelength to re
    bn = qgamma(0.5,2*n)
    re = h*(bn^n)
    return(re)
}

