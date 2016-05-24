## The Ecological Models and Data in R (EMD) Type-II curve
# Formally known as the Beddington-DeAngelis Type-II curve in FRAIR < 0.5
# Reduces to rogers_II when P == 1
# N0 replaced with 'X' for simplicity and consistency.
# X = Number of 'prey' (prey density / concentration)
# Y = Number of prey eaten / consumed / killed / absorbed

## EMD Type-II function ##
# Same as ?emdbook::lambertW, with the addition of 'P'
# Everything except 'X' should be provided.
emdII <- function(X, a, h, P, T) {
    if(is.list(a)){
        coefs <- a
        a <- coefs[['a']]
        h <- coefs[['h']]
        P <- coefs[['P']]
        T <- coefs[['T']]
    }
	return(X - lamW::lambertW0(a * h * X * exp(-a * (P * T - h * X)))/(a * h))
}

emdII_fit <- function(data, samp, start, fixed, boot=FALSE, windows=FALSE) {
	# Setup windows parallel processing
	fr_setpara(boot, windows)
	samp <- sort(samp)
	dat <- data[samp,]
	out <- fr_setupout(start, fixed, samp)

    try_emdII <- try(bbmle::mle2(emdII_nll, start=start, fixed=fixed, data=list('X'=dat$X, 'Y'=dat$Y), 
                         optimizer='optim', method='Nelder-Mead', control=list(maxit=5000)), 
                    silent=T)
	if (inherits(try_emdII, "try-error")) {
 		# The fit failed...
 		if(boot){
 			return(out)
        } else {
 			stop(try_emdII[1])
 		}
 	} else {
        # The fit 'worked'
 	    for (i in 1:length(names(start))){
            # Get coefs for fixed variables
 	        cname <- names(start)[i]
            vname <- paste(names(start)[i], 'var', sep='')
 	        out[cname] <- coef(try_emdII)[cname]
            out[vname] <- vcov(try_emdII)[cname, cname]
 	    }
 	    for (i in 1:length(names(fixed))){
            # Add fixed variables to the output
 	        cname <- names(fixed)[i]
 	        out[cname] <- as.numeric(fixed[cname])
 	    }
 		if(boot){
 			return(out)
 		} else {
 			return(list(out=out, fit=try_emdII))
 		}
 	}
}	
# emdII_nll
# Provides negative log-likelihood for estimations via bbmle::mle2()
# See Bowkers book for more info
emdII_nll <- function(a, h, P, T, X, Y) {
	if (a < 0 || h < 0){return(NA)}
	prop.exp = emdII(X, a, h, P, T)/X
	if(any(is.complex(prop.exp))){return(NA)} # Complex numbers don't help!
	# The proportion consumed must be between 0 and 1 and not NaN
	# If not then it must be bad estimate of a and h and should return NA
	if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)} 
	if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)} 
	return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}

# The difference function
emdII_diff <- function(X, grp, a, h, P, T, Da, Dh) {
  # return(X - lamW::lambertW0(a * h * X * exp(-a * (P * T - h * X)))/(a * h))
    return(X - lamW::lambertW0((a-Da*grp) * (h-Dh*grp) * X * exp(-(a-Da*grp) * (P * T - (h-Dh*grp) * X)))/((a-Da*grp) * (h-Dh*grp)))
}
# The diff NLL fucntions
emdII_nll_diff <- function(a, h, T, P, Da, Dh, X, Y, grp) {
    if (a < 0 || h < 0){return(NA)}
    prop.exp = emdII_diff(X, grp, a, h, P, T, Da, Dh)/X
    if(any(is.complex(prop.exp))){return(NA)} # Complex numbers don't help!
    # The proportion consumed must be between 0 and 1 and not NaN
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)} 
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)} 
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}



