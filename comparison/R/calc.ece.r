#############################################################################################
# functions to conduct emprical cross entropy calculations
# these functions are a reworking of Danial Ramos' ECE functions for Matlab
# (c) David Lucy 22 June 2010
# uses the gpava() function from the isotone() package as a calibrated
# set is needed
#############################################################################################
# calculates ECE according to Ramos (various)
# but mostly Ramos IEEE 2008 Equation 15
# has a default prior of 99 bins
calc.ece <- function(LR.ss, LR.ds, prior=seq(from=0.01, to=0.99, length=99))
{
# get the lengths of the LR vectors and prior vector
n.ss     <- length(LR.ss)
n.ds     <- length(LR.ds)
n.prior  <- length(prior)

# convert the prior to prior-odds
odds <- prior / (1 - prior)

# set up arrays of null LRs to calculate the 
# null ECE
LR.null.ss <- rep(1, n.ss)
LR.null.ds <- rep(1, n.ds)

# calculate a set of calibrated LRs
cal.set <- calibrate.set(LR.ss, LR.ds, method="raw")
LR.cal.ss <- cal.set$LR.cal.ss
LR.cal.ds <- cal.set$LR.cal.ds

# set up the empirical cross entropy vector
# and the null vector and calibrated array
ECE <- NULL
ECE.null <- NULL
ECE.cal <- NULL

	# for all values of the prior
	for(ctr in 1:n.prior)
		{
		bit.1 <- prior[ctr] / n.ss
		# for all same source LRs - do as a vector - sum later
		bit.2a <- log2(1 + (1 / (LR.ss * odds[ctr])))
		bit.2b <- log2(1 + (1 / (LR.null.ss * odds[ctr])))
		bit.2c <- log2(1 + (1 / (LR.cal.ss * odds[ctr])))

		bit.3 <- (1 - prior[ctr]) / n.ds
		# for all different source LRs - do as a vector - sum later
		bit.4a <- log2(1 + (LR.ds * odds[ctr]))
		bit.4b <- log2(1 + (LR.null.ds * odds[ctr]))
		bit.4c <- log2(1 + (LR.cal.ds * odds[ctr]))
	
		# the ECE for the ctr'th value of the prior
		# and the CE for the null for each value in the prior
		# Ramos IEEE 2008 Equation 15
		ECE[ctr] <- (bit.1 * sum(bit.2a)) + (bit.3 * sum(bit.4a))
		ECE.null[ctr] <- (bit.1 * sum(bit.2b)) + (bit.3 * sum(bit.4b))
		ECE.cal[ctr] <- (bit.1 * sum(bit.2c)) + (bit.3 * sum(bit.4c))
		}

# debug code
#assign("aa", prior.ratio[ctr], .GlobalEnv)

# S3 classes - function now S4
# put together a list of the prior and ECEs
# out <- list(prior, ECE, ECE.null, ECE.cal)
# names(out) <- c("prior", "ECE", "ECE.null", "ECE.cal")
# and return it
# return(out)

return(new("ece", prior=prior, ece.null=ECE.null, ece=ECE, ece.cal=ECE.cal))
}

