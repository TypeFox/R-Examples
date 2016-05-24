# Implements the elementary functions for a general
# step-up-down test. 
# Step-Up and step-down are derived from that.
# 
# Author: MarselScheer
###############################################################################


SUD <- function(pValues, criticalValues, startIDX_SUD) 
{		
	len <- length(criticalValues)
	
	# +++++++++++++++++   Plausis    ++++++++++++++
	if (len == 1)
		stop("SUD(): There is only 1 critical Value. Use the function SS()!")
	
	if (len != length(pValues))
		stop("SUD(): Length of critical Values and pValues need to be of the same length!")
	
	if (startIDX_SUD < 1 || len < startIDX_SUD )
		stop("SUD(): startIDX out of bound. criticalValues[startIDX] will not exist.")		
	# -----------------   Plausis    ---------------
	
	rejected <- rep(FALSE, times = len)

	# need to work with orderd pValues
	sortedPV <- sort(pValues, index.return = TRUE)
	
	suspiciousPV <- (sortedPV$x <= criticalValues)
	
	if (suspiciousPV[startIDX_SUD])
	{# Suspicious pValue.		
		# Actually doing now a StepDown on startIDX:len.
		# Additionally reject anything before startIDX
			
		nonSuspAboveStartIDX <- which(!suspiciousPV[startIDX_SUD:len]) 
		# ! looking only at the subset startIDX:len probably gives a shift of the index !
		# ! gonna correct this soon !
		
		# perhaps any pValue from startIDX to the end is suspicious, thus reject all!
		if (length(nonSuspAboveStartIDX) == 0)
			return(rep(TRUE, times = len))		
		
		# Correcting the shift
		nonSuspAboveStartIDX <- nonSuspAboveStartIDX + startIDX_SUD - 1
		
		# There must be some pValue between startIDX and the end that is not suspicious
		# Searching the first one. Anything immediately BEFORE that pValue will be rejected.		
		minIDX <- min(nonSuspAboveStartIDX) - 1
		
		rejected[sortedPV$ix[1:minIDX]] <- TRUE
	}
	else
	{# not suspicious pValue
		# Actually doing now a StepUp on 1:startIDX
		# The rejected are only the one rejected by this StepUp

		suspiciousIDX <- which(suspiciousPV[1:startIDX_SUD])
		
		# perhaps no pValue is suspicious, thus we do not reject anything
		if (length(suspiciousIDX) == 0)
			return(rep(FALSE, times = len))
		
		# There must be some pValue between 1 and startIDX that is suspicious
		# Searching the last one. Anything before (including the maximum) will be rejected.
		maxIDX <- max(suspiciousIDX)
		
		rejected[sortedPV$ix[1:maxIDX]] <- TRUE 
	}
	
	return(rejected)
}

SD <- function(pValues, criticalValues) 
{
	SUD(criticalValues = criticalValues, 
			pValues = pValues,
			startIDX_SUD = 1)
}

SU <- function(pValues, criticalValues) 
{
	SUD(criticalValues = criticalValues, 
			pValues = pValues,
			startIDX_SUD = length(criticalValues))
}



