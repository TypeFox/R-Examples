## Type I functional response function.
# A straight line with an intercept at zero, basically linear with respect to encouter rate

## Type I functional response ##
typeI <- function(X, a, T) {
    if(is.list(a)){
        coefs <- a
        a <- coefs[['a']]
        T <- coefs[['T']]
    }
    return(a*X*T) # Taken from Juliano 2001, eq. 10.1, pg 181. When h = 0 Ne=aNT/1
}

# typeI_fit: Does the heavy lifting
# data = The data from which to subsample. X and Y are drawn from here.
# samp = Provided by boot() or manually, as required
# start = List of starting values for items to be optimised.  Can only be 'a'.

typeI_fit <- function(data, samp, start, fixed, boot=FALSE, windows=FALSE) {
    # Setup windows parallel processing
	fr_setpara(boot, windows)
    
    samp <- sort(samp)
    dat <- data[samp,]
    out <- fr_setupout(start, fixed, samp)
    
    try_typeI <- try(bbmle::mle2(typeI_nll, start=start, fixed=fixed, data=list('X'=dat$X, 'Y'=dat$Y), optimizer='optim'), 
                     silent=T)
    ## Remove 'silent=T' for more verbose output
    if (inherits(try_typeI, "try-error")){
        # The fit failed...
        if(boot){
            return(out)
        } else {
            stop(try_typeI[1])
        }
    } else {
        # The fit 'worked'
        for (i in 1:length(names(start))){
            # Get coefs for fixed variables
            cname <- names(start)[i]
            vname <- paste(names(start)[i], 'var', sep='')
            out[cname] <- coef(try_typeI)[cname]
            out[vname] <- vcov(try_typeI)[cname, cname] 
        }
        for (i in 1:length(names(fixed))){
            # Add fixed variables to the output
            cname <- names(fixed)[i]
            out[cname] <- as.numeric(fixed[cname])
        }
        if(boot){
            return(out)
        } else {
            return(list(out=out, fit=try_typeI))
        }
    }
}	
# typeI_nll
# Provides negative log-likelihood for estimations via bbmle::mle2()
# See Bowkers book for more info
# Generalised from rogersII_nll, should be OK (DP)
typeI_nll <- function(a, T, X, Y) {
    if (a < 0) {return(NA)} # Zero would be a flat line, in this case, so is probably OK
    prop.exp = typeI(X, a, T)/X
    # The proportion consumed must be between 0 and 1 and not NaN
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)}  
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)} 
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}

# Type I difference function
typeI_diff <- function(X, grp, a, T, Da) {
  # return(a*X*T) # Taken from Juliano 2001, eq. 10.1, pg 181. When h = 0 Ne=aNT/1
    return((a-Da*grp)*X*T)
}

# The NLL for the difference model... used by frair_compare()
typeI_nll_diff <- function(a, T, Da, X, Y, grp) {
    if (a < 0) {return(NA)} # Zero would be a flat line, in this case, so is probably OK
    prop.exp = typeI_diff(X, grp, a, T, Da)/X
    # The proportion consumed must be between 0 and 1 and not NaN
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)}  
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)} 
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}
