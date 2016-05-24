BiCopTau2Par <- function(family, tau) {
    
    ## sanity checks
    if (length(family) != 1)
        stop("Input for family has to be a scalar/integer.")
    if (!(family %in% c(0, 1, 2, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36, 41, 51, 61, 71)))
        stop("Copula family not implemented.")
    
    
    ## calculation of parameter(s) depending on pair-copula family
    if (family == 0) {
        par <- rep(0, times = length(tau))
    } else if (family %in% 1:2) {
        par <- sin(pi * tau/2)
    } else if (family %in% c(3, 13)) {
        if (any(tau <= 0)) 
            stop("Clayton copula cannot be used for tau<=0.")
        par <- 2 * tau/(1 - tau)
    } else if (family %in% c(4, 14)) {
        if (any(tau < 0)) 
            stop("Gumbel copula cannot be used for tau<0.")
        par <- 1/(1 - tau)
    } else if (family == 5) {
        if (any(tau == 0)) 
            stop("Frank copula cannot be used for tau=0.")
        par <- sapply(tau, Frank.itau.JJ)
    } else if (family %in% c(6, 16)) {
        if (any(tau <= 0)) 
            stop("Joe copula cannot be used for tau<=0.")
        par <- sapply(tau, Joe.itau.JJ)
    } else if (family %in% c(23, 33)) {
        if (any(tau >= 0)) 
            stop("Rotated Clayton copula cannot be used for tau>=0.")
        par <- 2 * tau/(1 + tau)
    } else if (family %in% c(24, 34)) {
        if (any(tau > 0)) 
            stop("Rotated Gumbel copula cannot be used for tau>0.")
        par <- -(1/(1 + tau))
    } else if (family %in% c(26, 36)) {
        if (any(tau >= 0)) 
            stop("Rotated Joe copula cannot be used for tau>=0.")
        par <- -sapply(-tau, Joe.itau.JJ)
    } else if (family %in% c(41, 51)) {
        par <- sapply(tau, ipsA.tau2cpar)
    } else if (family %in% c(61, 71)) {
        par <- -sapply(-tau, ipsA.tau2cpar)
    }
    return(par)
    
}
