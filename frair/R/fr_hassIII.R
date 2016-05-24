## Hassell's Orginal Type III pre-prey function.
# hassIII: The guiding function...
hassIII <- function(X, b, c, h, T) {
    if(is.list(b)){
        coefs <- b
        b <- coefs[['b']]
        c <- coefs[['c']]
        h <- coefs[['h']]
        T <- coefs[['T']]
    }
    #a <- (d+b*X)/(1+c*X) # From Julliano (2001), pg. 181: a = (d+bN)/(1+cN)
    a <- (b*X)/(1+c*X) # From Hassell et al (1977)
    #cat(mean(a), b, c, h, '\n')
    return((a*X*T)/(1+a*X*h)) # Substituting into the code for hollingsII
	#return((d*X*T+b*X^2*T)/(1+c*X+d*X*h+b*X^2*h)) # Direct from Julliano (2001), pg 181
}
# hassIII_fit: Does the heavy lifting
hassIII_fit <- function(data, samp, start, fixed, boot=FALSE, windows=FALSE) {
	# Setup windows parallel processing
	fr_setpara(boot, windows)
	samp <- sort(samp)
	dat <- data[samp,]
	out <- fr_setupout(start, fixed, samp)

    try_hassIII <- try(bbmle::mle2(hassIII_nll, start=start, fixed=fixed, data=list('X'=dat$X, 'Y'=dat$Y), 
                            optimizer='optim', method='Nelder-Mead', control=list(maxit=5000)), 
                       silent=T)
	if (inherits(try_hassIII, "try-error")) {
 		# The fit failed...
 		if(boot){
 			return(out)
        } else {
 			stop(try_hassIII[1])
 		}
 	} else {
        # The fit 'worked'
 	    for (i in 1:length(names(start))){
            # Get coefs for fixed variables
 	        cname <- names(start)[i]
            vname <- paste(names(start)[i], 'var', sep='')
 	        out[cname] <- coef(try_hassIII)[cname]
            out[vname] <- vcov(try_hassIII)[cname, cname]
 	    }
 	    for (i in 1:length(names(fixed))){
            # Add fixed variables to the output
 	        cname <- names(fixed)[i]
 	        out[cname] <- as.numeric(fixed[cname])
 	    }
 		if(boot){
 			return(out)
 		} else {
 			return(list(out=out, fit=try_hassIII))
 		}
 	}
}	
# hassIII_nll: Provides negative log-likelihood for estimations via bbmle::mle2()
hassIII_nll <- function(b, c, h, T, X, Y) {
    if (h <= 0 || b <= 0) {return(NA)} # h and b estimates must be > zero
    if (c < 0) {return(NA)} # c must be positive (can be negative)
    prop.exp = hassIII(X, b, c, h, T)/X
    # The proportion consumed must be between 0 and 1 and not NaN
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)} 
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)} 
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}

# The diff function
hassIII_diff <- function(X, grp, b, c, h, T, Db, Dc, Dh) {
  # a <- (b*X)/(1+c*X) # From Hassell et al (1977)
    a <- ((b-Db*grp)*X)/(1+(c-Dc*grp)*X)
  # return((a*X*T)/(1+a*X*h)) # Substituting into the code for hollingsII
    return((a*X*T)/(1+a*X*(h-Dh*grp)))
    
}

# The diff_nll function
hassIII_nll_diff <- function(b, c, h, T, Db, Dc, Dh, X, Y, grp) {
    if (h <= 0 || b <= 0) {return(NA)} # h and b estimates must be > zero
    if (c < 0) {return(NA)} # c must be positive
    prop.exp = hassIII_diff(X, grp, b, c, h, T, Db, Dc, Dh)/X
    # The proportion consumed must be between 0 and 1 and not NaN
    # If not then it must be bad estimate of a and h and should return NA
    if(any(is.nan(prop.exp)) || any(is.na(prop.exp))){return(NA)} 
    if(any(prop.exp > 1) || any(prop.exp < 0)){return(NA)} 
    return(-sum(dbinom(Y, prob = prop.exp, size = X, log = TRUE)))
}