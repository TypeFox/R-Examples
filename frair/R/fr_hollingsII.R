## Holling's Orginal Type II pre-prey function.
# hollingsII: The guiding function...
hollingsII <- function(X, a, h, T) {
    if(is.list(a)){
        coefs <- a
        a <- coefs[['a']]
        h <- coefs[['h']]
        T <- coefs[['T']]
    }
	return((a*X*T)/(1+a*X*h)) # Direct from Julliano 2001, pp 181
}
# hollingsII_fit: Does the heavy lifting
hollingsII_fit <- function(data, samp, start, fixed, boot=FALSE, windows=FALSE) {
	# Setup windows parallel processing
	fr_setpara(boot, windows)
	samp <- sort(samp)
	dat <- data[samp,]
	out <- fr_setupout(start, fixed, samp)

    try_hollingsII <- try(bbmle::mle2(hollingsII_nll, start=start, fixed=fixed, data=list('X'=dat$X, 'Y'=dat$Y), 
                               optimizer='optim', method='Nelder-Mead', control=list(maxit=5000)), 
                          silent=T)
	if (inherits(try_hollingsII, "try-error")) {
 		# The fit failed...
 		if(boot){
 			return(out)
        } else {
 			stop(try_hollingsII[1])
 		}
 	} else {
        # The fit 'worked'
 	    for (i in 1:length(names(start))){
            # Get coefs for fixed variables
 	        cname <- names(start)[i]
            vname <- paste(names(start)[i], 'var', sep='')
 	        out[cname] <- coef(try_hollingsII)[cname]
            out[vname] <- vcov(try_hollingsII)[cname, cname]
 	    }
 	    for (i in 1:length(names(fixed))){
            # Add fixed variables to the output
 	        cname <- names(fixed)[i]
 	        out[cname] <- as.numeric(fixed[cname])
 	    }
 		if(boot){
 			return(out)
 		} else {
 			return(list(out=out, fit=try_hollingsII))
 		}
 	}
}	
# hollingsII_nll: Provides negative log-likelihood for estimations via bbmle::mle2()
hollingsII_nll <- function(a, h, T, X, Y) {
    if (a <= 0 || h <= 0) {return(NA)} # Estimates must be > zero
    prop.exp = hollingsII(X, a, h, T)/X
    # The proportion consumed must be between 0 and 1 and not NaN
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)} 
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)} 
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}

# The diff function
hollingsII_diff <- function(X, grp, a, h, T, Da, Dh) {
  # return(a*X*T/(1+a*X*h)) # Direct from Julliano 2001, pp 181
    return((a-Da*grp)*X*T/(1+(a-Da*grp)*X*(h-Dh*grp)))
}

# The diff_nll function
hollingsII_nll_diff <- function(a, h, T, Da, Dh, X, Y, grp) {
    if (a <= 0 || h <= 0) {return(NA)} # Estimates must be > zero
    prop.exp = hollingsII_diff(X, grp, a, h, T, Da, Dh)/X
    # The proportion consumed must be between 0 and 1 and not NaN
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)} 
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)} 
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}