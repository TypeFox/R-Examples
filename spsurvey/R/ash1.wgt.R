# File: ash1.wgt.r
# Purpose: Compute the average shifted histogram for weighted data
# Programmer: Tony Olsen
#	Based on original script by Susan Hornsby
# Date: February 1, 2005
# Revised: April 16, 2005
# Revised: July 29, 2005
# Revised: January 27, 2012

#
# ASH1 algorithm, to calculate the ASH estimate f[k] for k in 1:nbin.
#
ash1.wgt <- function(x, wgt=rep(1,length(x)), m=5, nbin=50, ab=NULL,
                     support="Continuous") {

# Bin the possibly weighted data
	v <- bin1.wgt(x, wgt, nbin, ab, support=support)

# Set delta based on range and number of bins	
	delta <- attr(v, "delta")
	h <- m*delta

# Set up vectors for estimation	
	nbin <- attr(v,"nbin")
	a <- attr(v,"ab")[1]
	b <- attr(v,"ab")[2]
  
# Add m-1 empty bins on ends when no ab boundary specified
	adj <- 0
	if( is.null(ab) ) { 
		v <- c(rep(0,m-1), v, rep(0,m-1))
		nbin <- nbin + 2* (m - 1)
		adj <- m - 1
	}
	else if( is.na(ab[1])) {
		v <- c(rep(0, m-1), v)
		nbin <- nbin + (m - 1)
		adj <- m-1
	}
	else if( is.na(ab[2])) {
		v <- c(v, rep(0,m-1) )
		nbin <- nbin
	}

# Compute lower limit, center, and upper limit of bins
	tlow <- a - adj*delta + ( (1:nbin) - 1.0)*delta
	tcen <- a - adj*delta + ( (1:nbin) - 0.5)*delta
	tup  <- a - adj*delta + (1:nbin)*delta

# Compute density height
	f <- rep(0,nbin)
	for(i in 1:nbin) {
		mlow <- max(1,i-m+1)
		mhi <- min(nbin,i+m-1)
		for(k in mlow:mhi) {
			if(tup[k] >= a & tlow[k] <= b ) 
				f[i] <- f[i] + v[k]*wgt.lim(k-i,m, mlow=mlow-i, mhi=mhi-i) 
		}
	}

# Adjust height so density area equals 1
	f <- f/(delta*sum(f))
  
# Construct output
	ash <- list(x=tcen, y=f)
	attr(ash, "delta") <- delta
	attr(ash, "m") <- m
	attr(ash, "h") <- h
	attr(ash, "support") <- support

# Return the result
	return(ash)
}


#
# BIN algorithm for unequal-probability sample, 
#
bin1.wgt <- function(x, wgt=rep(1,length(x)), nbin=50, 
                     ab=nicerange(x), support="Continuous") {

# Arguments
#   x  Vector of data to be used to estimate density. NAs are allowed.
#
#   wgt  vector of weights for each observation from a probability sample.
#            default is equal weights (equal probability)
#
#   nbin  Number of bins for density estimate
#
#   ab  optional range for support associated with the density.  Both values may
#       be equal to NA.  If equal to NA, then corresponding limit will be based
#       on nicerange().
#
#   support  If equal to "Continuous", then data are from a continuous
#            distribution.  If equal to "Ordinal", then data are from a discrete
#            distribution defined for integers only.

# Remove any missing data
	x <- x[!is.na(x)]
	wgt <- wgt[!is.na(x)]
	n <- length(x)

# Check that nbin is positive
	if (nbin <= 0) 
		stop("\nNumber of bin intervals nonpositive")

# Check for ab range
	tmp <- nicerange(x)
	if(is.null(ab)) {
		ab <- tmp
	} else {
		if( is.na(ab[1]) ) ab[1] <- tmp[1]
		if( is.na(ab[2]) ) ab[2] <- tmp[2]
		if (ab[1] >= ab[2])
			stop("\nInterval vector has negative orientation") 
	}

# Determine delta
# Continuous data case
	if(support == "Continuous") 
	delta <- (ab[2] - ab[1])/	nbin
	if(support == "Ordinal") {
		delta <- 1
		ab[1] <- floor(ab[1]) - 0.5
		ab[2] <- ceiling(ab[2]) + 0.5
		nbin <- ab[2] - ab[1] + 1
	}

# Sum weighted data in bins
	v <- rep(0,nbin)
	for(k in 1:nbin) {
		for(i in 1:n) {
			v[k] <- ifelse(((ab[1]+(k-1)*delta)<=x[i]) & (x[i]<(ab[1]+k*delta)), 
                     v[k] + wgt[i], v[k]) 
		}
	}
	attr(v, "nbin") <- nbin
	attr(v, "ab") <- ab
	attr(v, "delta") <- delta
	attr(v, "support") <- support

# Return the result
	return(v)
}


#
# Define weight function
#
wgt.lim <- function(i, m, mlow=(1-m), mhi=(m-1)) {
	K <- function(t) { (15/16)*(1-t^2)^2 }
	I <- mlow:mhi
	w <- m*K(i/m)/sum(K(I/m)) 
	return(w)
}


#
# Find nice range for binning
#
nicerange <- function (x, beta = 0.1) {
	ab <- range(x)
	del <- ((ab[2] - ab[1]) * beta)/2
	return(c(ab + c(-del, del)))
}
