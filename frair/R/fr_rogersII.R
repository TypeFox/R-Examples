## Roger's Type II decreasing prey function.
# Each function needs to have a specification (e.g. rogersII), a fit (e.g. rogersII_fit) and, where appropriate, a maximum likelihood NLL function.  
# Each function specification needs to be listed in 'resp_known' (in fr_functions.R) with a description.
# N0 replaced with 'X' for simplicity and consistency.
# X = Number of 'prey' (prey density / concentration)
# Y = Number of prey eaten / consumed / killed / absorbed

## Rogers Type II decreasing prey function ##
# Same as ?emdbook::lambertW
# Everything except 'X' should be provided.
rogersII <- function(X, a, h, T) {
    if(is.list(a)){
        coefs <- a
        a <- coefs[['a']]
        h <- coefs[['h']]
        T <- coefs[['T']]
    }
    return(X - lamW::lambertW0(a * h * X * exp(-a * (T - h * X)))/(a * h))

}
# rogersII_fit: Does the heavy lifting
# data = The data from which to subsample. X and Y are drawn from here.
# samp = Provided by boot() or manually, as required
# start = List of starting values for items to be optimised.  Usually 'a' and 'h'.
# fixed = List of 'Fixed data' (not optimised). Sometimes 'T', but I'm not too sure.
# Note required packages are reloaded here so Windows can do parallel computing!
# Not also that the statistic now (2013-04-13) now returns the variance 
rogersII_fit <- function(data, samp, start, fixed, boot=FALSE, windows=FALSE) {
	# Setup windows parallel processing
	fr_setpara(boot, windows)
	samp <- sort(samp)
	dat <- data[samp,]
	out <- fr_setupout(start, fixed, samp)

    try_rogersII <- try(bbmle::mle2(rogersII_nll, start=start, fixed=fixed, data=list('X'=dat$X, 'Y'=dat$Y), 
                             optimizer='optim', method='Nelder-Mead', control=list(maxit=5000)), 
                        silent=T) # Remove 'silent=T' for more verbose output
	if (inherits(try_rogersII, "try-error")) {
 		# The fit failed...
 		if(boot){
 			return(out)
        } else {
 			stop(try_rogersII[1])
 		}
 	} else {
        # The fit 'worked'
 	    for (i in 1:length(names(start))){
            # Get coefs for fixed variables
 	        cname <- names(start)[i]
            vname <- paste(names(start)[i], 'var', sep='')
 	        out[cname] <- coef(try_rogersII)[cname]
            out[vname] <- vcov(try_rogersII)[cname, cname]
 	    }
 	    for (i in 1:length(names(fixed))){
            # Add fixed variables to the output
 	        cname <- names(fixed)[i]
 	        out[cname] <- as.numeric(fixed[cname])
 	    }
 		if(boot){
 			return(out)
 		} else {
 			return(list(out=out, fit=try_rogersII))
 		}
 	}
}	
# rogersII_nll
# Provides negative log-likelihood for estimations via bbmle::mle2()
# See Ben Bowkers book for more info
rogersII_nll <- function(a, h, T, X, Y) {
    if (a <= 0 || h <= 0){return(NA)}
    prop.exp = rogersII(X, a, h, T)/X
    if(any(is.complex(prop.exp))){return(NA)} # Complex numbers don't help!
    # The proportion consumed must be between 0 and 1 and not NaN or NA
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)} 
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)}
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}

# Rogers II difference function
# Models the difference between two groups (j) exposing a simple t-test on Da and Dh
# For further info see Juliano 2001, pg 193, eg. eq. 10.11
rogersII_diff <- function(X, grp, a, h, T, Da, Dh) {
  # return(X - lamW::lambertW0(a * h * X * exp(-a * (T - h * X)))/(a * h))
    return(X - lamW::lambertW0((a-Da*grp) * (h-Dh*grp) * X * exp(-(a-Da*grp) * (T - (h-Dh*grp) * X)))/((a-Da*grp) * (h-Dh*grp))) 
}
# The NLL for the difference model... used by frair_compare()
rogersII_nll_diff <- function(a, h, T, Da, Dh, X, Y, grp) {
    if (a <= 0 || h <= 0){return(NA)}
    prop.exp = rogersII_diff(X, grp, a, h, T, Da, Dh)/X
    if(any(is.complex(prop.exp))){return(NA)} # Complex numbers don't help!
    # The proportion consumed must be between 0 and 1 and not NaN or NA
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)} 
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)}
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}
