# maRketSim - Core functions

# The heirarchy:
# core functions calculate values from values
# market objects calculates from core functions
# bond object and functions calculate from market object
# portfolio calculates from bond object and market object
# account calculates from bond object and market objects held within a portfolio

# 
# Time (in years):
# Each market object is a market in a single time period.  For convenience, we may want to create a market.history object which holds markets just like a portfolio holds bonds.
# Bonds inherit the market time when they are created and store that as their issue date.

# TODO
# Make accounts work - add history.account (just a data.frame of values over time) and plot.history.account
	# Current step: update.account needs to work
# Change "market" object class name to "market.bond"
# plot(history.market()) needs a legend
# Add zero-coupon bonds, stock markets, 
# Add daily pricing (between coupons)
# Add other kinds of risk premia to market objects (e.g. not just Treasuries)
# genPortfolio.bond periodically returns:
## Loop appears to have run away on us.
##Error in bond(mkt = mkt, mat = NA, dur = requested.dur, f = 0.5, ...) : 
##  (list) object cannot be coerced to type 'double'
#genHistory.market returns:
##Error in market(i = i, t = t, MM.frequency = f, ...) : 
##  unused argument(s) (i = i, MM.frequency = f)
# More efficient algorithm for max duration (1 over interest rate--confirm)
# Draw genPortfolio() from tnorm for more efficiency

# register package S3methods S3method(print,bond)


.onAttach <- function(libname = find.package("maRketSim"), pkgname = "maRketSim") {
	packageStartupMessage("maRketSim - This software is research software and is not intended for real-time use.  Under no circumstances should it be used as the basis for market or trading decisions.")
}




# --- Core functions --- #

# Find maturity from duration
# Right now this rounds to the nearest coupon date
findMat <- function(dur,i,f=.5) {
	min.search <- ifelse(floor(dur)>0,floor(dur),.5) # ensure we can't be starting at a maturity of zero
  max.search <- dur*10 # a better solution, but one which leads to circularities since it depends on this function, is findMaxMat(i=i,market.rate=i)
	ms <- seq(min.search,max.search,f)
	ds <- sapply(ms,findDur,i=i,f=f)
  ms.sel <- which.min(abs(ds-dur))
  if( ms.sel==length(ms) ) { warning("findMat hit max.search, which likely means there was a problem.\n Maybe you gave the function a duration greater than max duration under that interest rate?\n") }
	return(ms[ms.sel])
}

# Find duration from maturity
findDur <- function(mat,i,market.rate=NA,f=.5,type="modified",...) {
	# Error check inputs
	if(type!="modified"&type!="Macaulay") stop("Must specify modified or Macaulay for record type.\n")
	# Set numbers which cancel out
	par=1000
	cpn=par*i*f
	# Calculate intermediate values
	if(is.na(market.rate)) market.rate=i
	pv <- findPV(i=i,market.rate=market.rate,mat=mat,f=f,par=par)
	# Calculate periods based on whether or not the bond is on its coupon date
	if(((mat/f)-as.integer(mat/f))!=0) { # If the bond isn't on its coupon date
		ms <- seq(as.integer(mat),as.integer(mat)+1,f)
		mat.whole <- ms[(mat-ms)>=0][length(ms[(mat-ms)>=0])]
		mat.frac <- mat-mat.whole
		w <- mat.frac/f # fraction of the period remaining until the next coupon payment
		n <- mat.whole/f # number of whole coupon payments remaining
		periods <- seq(1,n)-1+w
	} else { # If the bond falls on its coupon date
		n <- mat/f
		periods <- seq(1,n)
	}
	# Duration calculation
	dur.Macaulay <- f * (sum(periods * pv$pv.cf) / pv$pv ) # Macaulay duration, in years
	if(type=="Macaulay") {
		cat("Returning Macaulay duration\n")
		return(dur.Macaulay)
	} else {
		dur <- dur.Macaulay / (1+market.rate*f)
		return(dur)
	}
}

# Find PV of a bond and its cash flows - the fundamental function used by many others (including constructing a bond object)
findPV <- function(i,mat,market.rate,par=1000,f=.5,fractional.method="30/360") {
	if(((mat/f)-as.integer(mat/f))!=0) stop("No between-coupon maturities please.  This is currently producing discontinuous values in an odd way.")
	# Error check inputs
	if(fractional.method!="30/360") stop("actual/actual mid-period accounting not yet supported") # see Handbook of Fixed Income, 1991, p91
	# Calculate intermediate values
	cpn <- par*i*f
	# When settlement date falls between coupon periods
	if(fractional.method=="30/360" & ((mat/f)-as.integer(mat/f))!=0) {
		# Find fractional part of the maturity
		ms <- seq(as.integer(mat),as.integer(mat)+1,f)
		mat.whole <- ms[(mat-ms)>=0][length(ms[(mat-ms)>=0])]
		mat.frac <- mat-mat.whole
		w <- mat.frac/f # fraction of the period remaining until the next coupon payment
		n <- mat.whole/f # number of whole coupon payments remaining
		ws <- c(w,rep(1,n-1)) # Multiply only the first term by w
		periods <- w+seq(1,n)-1
		# Calculate present value of each cash flow
		pv.cpn <- ws*cpn/(1+market.rate*f)^periods
		pv.par <- par/(1+market.rate*f)^(n-1+w)
	} else { # When settlement date falls on a coupon period
		n <- mat/f # number of periods
		periods <- seq(1,n)
		pv.cpn <- cpn/(1+market.rate*f)^periods
		pv.par <-   par/(1+market.rate*f)^periods[n]
	}
	# - Calculate PV - #
	pv.cf <- pv.cpn
	pv.cf[n] <- pv.cpn[n] + pv.par
	pv <- sum(pv.cf)
	# Return results
	res <- list(pv=pv, pv.cpn = pv.cpn, pv.par = pv.par, pv.cf = pv.cf)
	res
}

# Find the future value of a cashflow
findFV <- function(P,i,t.elapsed,compound="continuous") {
	if(compound=="continuous") {		
		fv <- P * exp(i*t.elapsed)
	} else if (!is.na(as.numeric(compound))) {
		n <- 1/as.numeric(compound)
		fv <- P * (1+i/n)^(t.elapsed*n)
	} else { stop("Must specify compound as continuous or as a frequency e.g. 0.5 is semi-annual compounding.\n") }
	fv
}

# Estimate duration using various closed-form formulae
# Equations 5,6, and 11 in Pianca "Maximum duration of below-par bonds: A closed-form formula"
# Assumptions: flat yield curve, constant coupon, reimbursement value = face value
# r = C/F, or the coupon rate (F is face value, C is dollar value of coupons)
# i = applicable interest rate
# n = maturity date,
# type = "pianca", "macaulay", or "hawawini"
# At par, (i==r)
findDur_ClosedForm <- function(i,market.rate,mat,type="pianca") {
  type <- tolower(type)
  ani <- NA # For hawawini: Need pv of an n-period annuity at rate i
  switch(type,
    pianca = 1 + (1/market.rate) + ( mat*(market.rate-i) - (1+market.rate) )/( i*( (1+market.rate)^mat - 1 ) + market.rate ) ,
    macaulay = 1 + 1/market.rate - ( (1+market.rate)/i + mat*(1+1/i-(1+market.rate)/i ) ) / ( (1+market.rate)^mat - 1 - 1/i + (1+market.rate)/i ) ,
    hawawini = ( (1+market.rate)*ani*i + mat*(market.rate-i)(1+market.rate)^(-mat) ) / ( i+(market.rate-i)(1+market.rate)^-mat )
  )
}

# Find maximum duration using closed-form formulae
# ... pass-alongs to findDur_ClosedForm
# i is the coupon rate
# market.rate is the current market rate
findMaxDur <- function( i, market.rate, ... ) {
  # If above or at par, max duration is 1+1/i
  # Otherwise use Pianca formula
  asymptote <- 1+1/market.rate
  if( market.rate <= i ) { # At or above par
    return(asymptote)
  } else { # Below par
    n_star <- findMaxMat( i=i, market.rate=market.rate, ... )
    return( findDur_ClosedForm( i=i, market.rate=market.rate, mat=n_star ) )
  }
}

findMaxMat <- function( i, market.rate, ... ) {
  if( market.rate <= i ) { # At or above par
    asymptote <- 1+1/market.rate
    return( findMat(dur=asymptote,i=i) )
  } else { # Below par
    require(gsl) # load Lambert's W function
    a <- market.rate-i
    b <- log(1+market.rate)
    n_star <- ( b*(1+market.rate) + a*( 1 + lambert_W0( a*exp( -(a+b*(1+market.rate))/a )/i ) ) ) / (a*b)
    return( n_star )
  }
}
