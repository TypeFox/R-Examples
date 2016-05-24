convrad = function(n, f = 0.9, r = 1, fr = 0.5){
    # converts sersic radii containing differing amounts of total luminosity
    bn1 = qgamma(fr,2*n)
    bn2 = qgamma(f,2*n)
    r2 = r*((bn2/bn1)^n)
    return(r2)
}

