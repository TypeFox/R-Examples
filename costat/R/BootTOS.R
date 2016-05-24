BootTOS <-
function(x, Bsims=100, WPsmooth=TRUE, verbose=FALSE, plot.avspec=FALSE,
	plot.avsim=FALSE, theTS = TOSts, AutoReflect=TRUE, lapplyfn=lapply){
#
# Compute bootstrap test of stationarity on time series x
#
#
# Get the name of the time series
#
DNAME <- deparse(substitute(x))
#
# For time series, one way of dealing with boundaries is reflection
#
if (AutoReflect==TRUE)	{
	x <- c(x, rev(x))
	}

#
# Compute EWS estimate of x
#
xews <- ewspec(x, WPsmooth=WPsmooth, smooth.dev=var)$S

if (AutoReflect==TRUE)	{
	xews <- AntiAR(xews)
	x <- x[1:(length(x)/2)]
	}
#
# Compute Test Statistic for the data, x
#
TSdata <- theTS(xews)
#
# Now compute mean spectrum (or the null spectrum)
#
J <- xews$nlevels
n <- length(x)
m <- matrix(xews$D, nrow=J, ncol=n, byrow=TRUE) # Turn spec into matrix
m <- apply(m, 1, mean) # Vector of means one for each scale
m[m < 0] <- 0
m <- matrix(m, nrow=J, ncol=n)
xavspec <- xews
xavspec$D <- as.vector(t(m))

if (plot.avspec==TRUE)
	plot(xavspec)
#
# Now do bootstrap simulations. Generate time series with spectrum xavspec,
# compute test statistic and store.
#
# The following internal function simulates from the average spectrum in
# xavspace, computes its EWS, applies the test statistic and returns the
# result.
#

bsfn <- function(dummy, xavspec, WPsmooth, smooth.dev, theTS){

	xbs <- LSWsim(xavspec)
	xbs.ews <- ewspec(xbs, WPsmooth=WPsmooth, smooth.dev=smooth.dev)$S
	ans <- theTS(xbs.ews)
	return(ans)
}
dummy.ip <- vector("list", Bsims-1)

ans <- lapplyfn(dummy.ip, bsfn, xavspec=xavspec, WPsmooth=WPsmooth,
	smooth.dev=var, theTS=theTS)

TS <- unlist(ans)

TS <- c(TSdata, TS)

if (verbose==TRUE)
	cat("\n")

#
# Compute the p-value of the test
#
p.value <- plotBS(TS, plot=FALSE, verbose=FALSE)

#
# Return the results of this hypothesis test as an htest object
#
htest.obj <- list(statistic=TS[1], p.value=p.value, method="BootTOS test of stationarity", data.name= DNAME, Bootvals=TS)
#
# The returned object should inherit from htest and be of primary class
# BootTOS
#
class(htest.obj) <- c("BootTOS", "htest")

return(htest.obj)
}
