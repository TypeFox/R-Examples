#############################################################################################
# functions to conduct emprical cross entropy calculations
# these functions are a reworking of Danial Ramos' ECE functions for Matlab
# (c) David Lucy 20 July 2010
# calculates and returns the calibrated set of `ideal' LRs from the observed
# LRs using the penalised adjacent violators algotrthm
# This is very much a rewrite of Nico Brummer's optloglr() function
# for Matlab, only for R
# draws heavily on Nico Brummer's work
#############################################################################################
calibrate.set <- function(LR.ss, LR.ds, method="raw")
{
# isotone() needed for gpava()
# require(isotone) # require not needed as isotone already loaded

# correction for sorting - possibly not needed for R's sorting
# leave it in anyway
LR.ss <- LR.ss - 1E-6

# get the lengths of the LR vectors and prior vector
n.ss     <- length(LR.ss)
n.ds     <- length(LR.ds)
n  	 <- n.ss + n.ds

# set an arbitary cutoff to prevent division by zero errors
# but bodgy this and I'm not too keen on it - however
# so long as the value is close to, but less than, 1 it makes little difference
# cutoff <- 0.99999999999999999999999999

# produce vectors for the indicator function and LRs
# use convention that 1 signifies same set, and 0 different set
indicator.array <- c(rep(0, length=n.ds), rep(1, length=n.ss))
LR.array <- c(LR.ds, LR.ss)

# order the vectors
ordered.indicies <- order(LR.array)
ordered.indicator.array <- indicator.array[ordered.indicies]

# implement Laplace's rule of sucession in the way of Brummer
if(method == "laplace"){ordered.indicator.array <- c(1, 0, ordered.indicator.array, 1, 0)}


# gpava doesn't really care what you send it as the first arguement so long as
# it is ascending and of the appropriate length - here I just sent an integer array
if(method == "raw"){calibrated.set <- gpava(1:n, ordered.indicator.array)}
if(method == "laplace"){calibrated.set <- gpava(1:(n+4), ordered.indicator.array)}

calibrated.posterior.probabilities <- calibrated.set$x

# disentangle all the bits and pieces
if(method == "laplace")
	{
	calibrated.posterior.probabilities <- calibrated.posterior.probabilities[3:(n+2)]
	ordered.indicator.array <- ordered.indicator.array[3:(n+2)]
	}


# prior odds 
prior.odds <- n.ss / n.ds

# arbitrarily set a cutoff to prevent division by zero errors when
# we calculate the odds not very keen on this - don't think it is needed
# calibrated.posterior.probabilities[calibrated.posterior.probabilities > cutoff] <- cutoff

log.calibrated.posterior.LRs <- log(calibrated.posterior.probabilities / (1 - calibrated.posterior.probabilities)) - log(prior.odds)

# Brummer adds this in the ensure the idempotent property of the logLRs - not too sure how important this is
bit.to.add.on <- 1:n * 1E-6 / n
log.calibrated.posterior.LRs <- log.calibrated.posterior.LRs + bit.to.add.on
calibrated.posterior.LRs <- exp(log.calibrated.posterior.LRs)

# unpack all the calibrated LR values
LR.cal.ss <- calibrated.posterior.LRs[ordered.indicator.array == 1]
LR.cal.ds <- calibrated.posterior.LRs[ordered.indicator.array == 0]

out <- list(LR.cal.ss, LR.cal.ds)
names(out) <- c("LR.cal.ss", "LR.cal.ds")

# debug code
# assign("kk", calibrated.posterior.LRs, .GlobalEnv)

return(out)
}

