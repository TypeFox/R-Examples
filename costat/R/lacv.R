lacv <-
function(x, filter.number=10, family=c("DaubExPhase","DaubLeAsymm"),
	smooth.dev=var, AutoReflect=TRUE, lag.max=NULL, smooth.RM=0, ...)
{
family <- match.arg(family)

#
# Computes localized autocovariance/autocorrelation of time series x
#

#
# Better results are often achieved with AutoReflection as it constructs
# an artificial time series with periodic boundaries which matches the
# transforms used.
#

if (AutoReflect==TRUE)	{
	x <- c(x, rev(x))
	}

#
# Choose smoothing method
#

if (smooth.RM == 0)	# Wavelet smoothing a la NvSK00
	xews <- ewspec(x, filter.number=filter.number, family=family,  smooth.dev=smooth.dev, ...)$S
else	{
	#
	# Linear running mean smoothing with bandwidth smooth.RM
	#
	xews <- ewspec(x, filter.number=filter.number, family=family,  WPsmooth=FALSE, ...)$S
	xews <- EWSsmoothRM(xews, s=smooth.RM)
	}

#
# Get dimensions and turn spectrum into a matrix object
#

J <- xews$nlevels
xewsm <- matrix(xews$D, nrow=length(x), ncol=J)
#
# First col of xewsm is finest
#
# Rows are time
#
# If we did AutoReflect then we can get rid of second half of rows
# Also get rid of coarsest column (ie far right)
#
if (AutoReflect==TRUE)	{
	nr <- nrow(xewsm)
	xewsm <- xewsm[1:(nr/2), 1:(ncol(xewsm)-1)]
	J <- J - 1
	}

#
# Compute (nondecimated) autocorrelation wavelets using wavethresh,
# and related dimensions
#

Psi <- PsiJmat(-J, filter.number=filter.number, family=family)
nc <- ncol(Psi)
L <- (nc-1)/2
dimnames(Psi) <- list(NULL, c(-L:0, 1:L))

#
# Only need columns 0:L, as it's symmetric
#
if (is.null(lag.max))	# Use all lags
	the.lacv <- xewsm %*% Psi[, (L+1):ncol(Psi)]
else	{
	#
	# Use bespoke number of lags
	#
	if (L+1+lag.max > ncol(Psi)) 	{
		warning(paste("lag.max too high. Have reset it to ",
			ncol(Psi)-L-1, ". Higher lags are zero"))
		lag.max <- ncol(Psi)-L-1
		}
	the.lacv <- xewsm %*% Psi[, (L+1):(L+1+lag.max)]
	}

#
# Compute autocorrelations efficiently using autocovariance using sweep
#

the.lacor <- sweep(the.lacv, 1, the.lacv[,1], FUN="/")

#
# Construct return object, make it the right class and return it
#

l <- list(lacv=the.lacv, lacr=the.lacor, date=date())
class(l) <- "lacv"
return(l)
}
