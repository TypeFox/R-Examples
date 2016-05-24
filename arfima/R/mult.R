"shift" <- function(pars, by.am, useCt = TRUE) {
    
    if (useCt && by.am > 0) {
        ss <- function(vv, byam) .C("shift_C", as.double(vv), as.integer(byam), as.integer(length(vv)), 
            y = double((length(vv) - 1) * byam + 1))$y
        return(ss(pars, by.am))
    }
    n <- length(pars)
    
    if (all(pars[-1] == 0)) 
        return(pars)
    
    if (by.am < 0) {
        by.am <- abs(by.am)
        m <- floor((n - 1)/by.am)
        
        seq.pol <- rep(0, m + 1)
        
        seq.pol[1] <- pars[1]
        
        for (i in 2:(m + 1)) seq.pol[i] <- pars[(i - 1) * by.am + 1]
        
        return(seq.pol)
        
    }
    
    seq.pol <- rep(0, (n - 1) * by.am + 1)
    
    seq.pol[1] <- pars[1]
    
    for (i in 2:n) seq.pol[(i - 1) * by.am + 1] <- pars[i]
    
    return(seq.pol)
}

"mult" <- function(a, b) {
    if (all(a == 0) || all(b == 0)) 
        return(0)
    ret <- convolve(a, rev(b), type = "open")
    return(ret)
}

expand = function(d, Bterm = -1, n = 10) {
    
    ret = as.vector(sapply(0:n, function(x) choose(d, x) * Bterm^x))
    
    return(ret)
    
}

expandseas = function(d, Bterm = -1, seas = 12, n = 10) {
    
    return(shift(expand(d, Bterm, n), seas))
    
} 
