## The scaling exponent formualtion originally(?) stemming from uses of Hall's exponents by Real (1977)
# flexp: A typeII + scaling exponent
flexp <- function(X, b, q, h, T) {
    if(is.list(b)){
        coefs <- b
        b <- coefs[['b']]
        q <- coefs[['q']]
        h <- coefs[['h']]
        T <- coefs[['T']]
    }
    a <- (b*X^q)
    return((a*X*T)/(1+a*X*h))
    ## Original / real77 (incorrect)
    #a <- (b*X^q)
    #return((a*X^(q+1)*T)/(1+a*X^(q+1)*h))
}

# flexp_fit: Does the heavy lifting
flexp_fit <- function(data, samp, start, fixed, boot=FALSE, windows=FALSE) {
    # Setup windows parallel processing
    fr_setpara(boot, windows)
    samp <- sort(samp)
    dat <- data[samp,]
    out <- fr_setupout(start, fixed, samp)
    
    try_flexp <- try(bbmle::mle2(flexp_nll, start=start, fixed=fixed, data=list('X'=dat$X, 'Y'=dat$Y), 
                           optimizer='optim', method='Nelder-Mead', control=list(maxit=5000)), 
                      silent=T)
    if (inherits(try_flexp, "try-error")) {
        # The fit failed...
        if(boot){
            return(out)
        } else {
            stop(try_flexp[1])
        }
    } else {
        # The fit 'worked'
        for (i in 1:length(names(start))){
            # Get coefs for fixed variables
            cname <- names(start)[i]
            vname <- paste(names(start)[i], 'var', sep='')
            out[cname] <- coef(try_flexp)[cname]
            out[vname] <- vcov(try_flexp)[cname, cname]
        }
        for (i in 1:length(names(fixed))){
            # Add fixed variables to the output
            cname <- names(fixed)[i]
            out[cname] <- as.numeric(fixed[cname])
        }
        if(boot){
            return(out)
        } else {
            return(list(out=out, fit=try_flexp))
        }
    }
}

# flexp_nll: Provides negative log-likelihood for estimations via bbmle::mle2()
flexp_nll <- function(b, q, h, T, X, Y) {
    if (b <= 0 || h <= 0) {return(NA)} # Estimates must be > zero
    if (q < -1){return(NA)} # q+1 must be positive
    prop.exp = flexp(X, b, q, h, T)/X
    # The proportion consumed must be between 0 and 1 and not NaN
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)} 
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)} 
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}

# The diff function
flexp_diff <- function(X, grp, b, q, h, T, Db, Dq, Dh) {
  # a <- ( b        *X^ q        )
    a <- ((b-Db*grp)*X^(q-Dq*grp))
  # return((a*X*T)/(1+a*X* h        ))
    return((a*X*T)/(1+a*X*(h-Dh*grp)))
}

# The diff_nll function
flexp_nll_diff <- function(b, q, h, T, Db, Dq, Dh, X, Y, grp) {
    if (b <= 0 || h <= 0) {return(NA)} # Estimates must be > zero
    if (q < -1){return(NA)} # q+1 must be positive
    prop.exp = flexp_diff(X, grp, b, q, h, T, Db, Dq, Dh)/X
    # The proportion consumed must be between 0 and 1 and not NaN
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)} 
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)} 
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}

