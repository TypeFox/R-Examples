## Hassell's Type III pred-prey function not assuming replacement.
# hassIIInr: The guiding function...
hassIIInr <- function(X, b, c, h, T) {
    if(is.list(b)){
        coefs <- b
        b <- coefs[['b']]
        c <- coefs[['c']]
        h <- coefs[['h']]
        T <- coefs[['T']]
    }
    a <- (b*X)/(1+c*X) # From Hassell et al (1977)
    return(X - lamW::lambertW0(a * h * X * exp(-a * (T - h * X)))/(a * h)) # Direct from rogersII

}
# hassIIInr_fit: Does the heavy lifting
hassIIInr_fit <- function(data, samp, start, fixed, boot=FALSE, windows=FALSE) {
	# Setup windows parallel processing
	fr_setpara(boot, windows)
	samp <- sort(samp)
	dat <- data[samp,]
	out <- fr_setupout(start, fixed, samp)

	try_hassIIInr <- try(bbmle::mle2(hassIIInr_nll, start=start, fixed=fixed, data=list('X'=dat$X, 'Y'=dat$Y), 
	                         optimizer='optim', method='Nelder-Mead', control=list(maxit=5000)), 
	                   silent=T)
	if (inherits(try_hassIIInr, "try-error")) {
 		# The fit failed...
 		if(boot){
 			return(out)
        } else {
 			stop(try_hassIIInr[1])
 		}
 	} else {
        # The fit 'worked'
 	    for (i in 1:length(names(start))){
            # Get coefs for fixed variables
 	        cname <- names(start)[i]
            vname <- paste(names(start)[i], 'var', sep='')
 	        out[cname] <- coef(try_hassIIInr)[cname]
            out[vname] <- vcov(try_hassIIInr)[cname, cname]
 	    }
 	    for (i in 1:length(names(fixed))){
            # Add fixed variables to the output
 	        cname <- names(fixed)[i]
 	        out[cname] <- as.numeric(fixed[cname])
 	    }
 		if(boot){
 			return(out)
 		} else {
 			return(list(out=out, fit=try_hassIIInr))
 		}
 	}
}	
# hassIIInr_nll
hassIIInr_nll <- function(b, c, h, T, X, Y) {
    if (h <= 0 || b <= 0) {return(NA)} # h and b estimates must be > zero
    if (c < 0) {return(NA)} # c must be positive
    prop.exp = hassIIInr(X, b, c, h, T)/X
    if(any(is.complex(prop.exp))){return(NA)} # Complex numbers don't help!
    # The proportion consumed must be between 0 and 1 and not NaN
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)} 
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)} 
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}

# hassIIInr_diff
# Model the difference between two groups (j) exposing a simple t-test on Da and Dh
# For further info see Juliano 2001, pg 193, eg. eq. 10.11
hassIIInr_diff <- function(X, grp, b, c, h, T, Db, Dc, Dh) {
  # a <- ( b        *X)/(1+ c        *X) # From Hassell et al (1977)
    a <- ((b-Db*grp)*X)/(1+(c-Dc*grp)*X)
  # return(X-lamW::lambertW0(a* h        *X*exp(-a*(T-h         *X)))/(a* h))
    return(X-lamW::lambertW0(a*(h-Dh*grp)*X*exp(-a*(T-(h-Dh*grp)*X)))/(a*(h-Dh*grp))) 
}

# The NLL for the difference model... used by frair_compare()
hassIIInr_nll_diff <- function(b, c, h, T, Db, Dc, Dh, X, Y, grp) {
    if (h <= 0 || b <= 0) {return(NA)} # h and b estimates must be > zero
    if (c < 0) {return(NA)} # c must be positive
    prop.exp = hassIIInr_diff(X, grp, b, c, h, T, Db, Dc, Dh)/X
    if(any(is.complex(prop.exp))){return(NA)} # Complex numbers don't help!
    # The proportion consumed must be between 0 and 1 and not NaN
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)} 
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)} 
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}
