"psiwtsAR" <- function(phi, maxlag) {
    p <- length(phi)
    if (p == 0) 
        return(0)
    x <- numeric(maxlag + 1)
    x <- 1
    for (i in 1:p) {
        x[i + 1] <- crossprod(phi[1:i], (rev(x))[1:i])
    }
    if (maxlag > p) {
        for (i in (p + 1):maxlag) {
            x[i + 1] <- crossprod(phi, (rev(x))[1:p])
        }
    }
    return(x)
}


## Not used right now.
piwts <- function(phi = 0, theta = 0, phiseas = 0, thetaseas = 0, dfrac = 0, dfs = 0, dint = 0, 
    dseas = 0, period = 0, n = 128, m = 128, maxlag = 128) {
    multer = 1
    
    if ((length(phi) > 0) && any(phi != 0)) 
        multer <- mult(multer, c(1, -phi))
    
    if ((length(phiseas) > 0) && any(phi != 0) && (period > 0)) 
        multer <- mult(multer, shift(c(1, -phiseas), period))
    
    if ((length(dfrac) > 0) && dfrac != 0) 
        multer <- mult(multer, expand(d = dfrac, n = n))
    
    if ((length(dfs) > 0) && dfs != 0) 
        multer <- mult(multer, expandseas(d = dfs, seas = period, n = n * period))
    
    if ((length(dint) > 0) && dint > 0) {
        multerd <- 1
        for (i in 1:dint) multerd <- mult(multerd, c(1, -1))
        multer <- mult(multer, multerd)
    }
    
    if ((length(dint) > 0) && (dseas > 0) && (period > 0)) {
        multerds <- 1
        for (i in 1:dseas) multerds <- mult(multerds, shift(c(1, -1), period))
        multer <- mult(multer, multerds)
    }
    
    if ((length(theta) > 0) && any(theta != 0)) 
        multer <- mult(multer, psiwtsAR(theta, m))
    
    if ((length(thetaseas) > 0) && any(thetaseas != 0) && (period > 0)) 
        multer <- mult(multer, psiwtsAR(shift(c(1, thetaseas), period)[-1], m * period))
    
    return(multer[1:maxlag])
}


psiwts <- function(phi = 0, theta = 0, phiseas = 0, thetaseas = 0, dfrac = 0, dfs = 0, dint = 0, 
    dseas = 0, period = 0, len = 128, n = len * 2, div = 1) {
    multer <- 1
    if ((length(theta) > 0) && any(theta != 0)) 
        multer <- mult(multer, c(1, -theta))
    if ((length(phi) > 0) && any(phi != 0)) 
        multer <- mult(multer, psiwtsAR(phi, n))
    if (length(dint) > 0 && dint > 0 && length(dfrac) > 0 && dfrac != 0) 
        multer <- mult(multer, expand(d = -(dint + dfrac), n = n)) else if (length(dint) > 0 && dint > 0) 
        multer <- mult(multer, expand(d = -dint, n = n)) else if (length(dfrac) > 0 && dfrac != 0) 
        multer <- mult(multer, expand(d = -dfrac, n = n))
    
    if (period > 0) {
        if (length(phiseas) > 0 && any(phiseas != 0)) 
            multer <- mult(multer, psiwtsAR(shift(c(1, phiseas), period)[-1], n * period))
        if (length(thetaseas) > 0 && any(thetaseas != 0)) 
            multer <- mult(multer, shift(c(1, -thetaseas), period))
        if (length(dseas) > 0 && dseas > 0 && length(dfs) > 0 && dfs != 0) 
            multer <- mult(multer, expandseas(d = -(dseas + dfs), seas = period, n = n)) else if (length(dseas) > 0 && dseas > 0) 
            multer <- mult(multer, expandseas(d = -dseas, seas = period, n = ceiling(n/div))) else if (length(dfs) > 0 && dfs != 0) 
            multer <- mult(multer, expandseas(d = -dfs, seas = period, n = n))
    }
    if (length(multer) < len) 
        multer <- c(multer, rep(0, len))
    return(multer[1:len])
}

wtsforexact <- function(dint = 0, dseas = 0, period = 0, len) {
    if (dint != round(dint) || dseas != round(dseas) || period != round(period)) 
        stop("non-integers provided where integers needed")
    if (length(dseas) > 0 && dseas > 0 && period < 2) 
        stop("incompatible period and seasonal d")
    
    if (length(dint) > 0 && dint > 0) 
        nonseas <- sapply(0:(len - 1), function(j) (-1)^j * choose(-dint, j))
    
    if (period > 0 && length(dseas) > 0 && dseas > 0) {
        seas <- shift(sapply(0:(len - 1), function(j) (-1)^j * choose(-dseas, j)), period)
    } else seas <- 1
    mult <- mult(nonseas, seas)  ##this is what makes this approximate.
    return(mult[1:len])
} 
