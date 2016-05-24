################################################################################
# Two functions for checking and fixing tied event-times. 

################################################################################
# Check whether there are tied event-times, and if there are and 
# precision!=NULL then deal with them. Example: precision=1/365.25 means times
# were recorded to the nearest day but are stored as years. If precision!=NULL 
# then also correct the times, by rounding them to the nearest day or whatever.  

fixTiedEventTimes <- function(timeAtEntry, timeAtExit, isCase, precision) {
	
	# If precision is not null, and there are tied event-times, then fix them.
	if (!is.null(precision)) {
		
		# Check that timeAtEntry/Exit are recorded to the nearest "precision". 
		shouldBeIntegers <- c(timeAtEntry, timeAtExit) / precision
		if (any(abs(shouldBeIntegers - round(shouldBeIntegers)) > 0.01)) {
			if (precision==1/365.25) precision <- "1/365.25"
			if (precision==1/12) precision <- "1/12"
			if (precision==1/7) precision <- "1/7"
			stop("some entry/exit times do not seem to have been ",
					"measured to the nearest precision=", precision)
		}
		
		# Create integer versions of the times.
		timeAtEntryInteger <- as.integer(round(timeAtEntry/precision))
		timeAtExitInteger <- as.integer(round(timeAtExit/precision))
		
		# Change any tied integer-event-times to get rid of ties:
		perturbResult <- perturbTiedValues(
				timeAtExitInteger[isCase], epsilonForTies=0.01)
		timeAtExitInteger[isCase] <- perturbResult$values
		message <- perturbResult$message
		# (What should epsilonForTies be? The values are integers, so it 
		# certainly needs to be less than 1. In non-thorough experiments, 0.01
		# and 0.001 lead to the same values in the eventual output from cchs.) 
		
		# Replace original times by accurate and corrected times:
		timeAtEntry <- timeAtEntryInteger * precision
		timeAtExit <- timeAtExitInteger * precision
		
	} else {
		message <- ""
	}
	
	# Now check whether there are any duplicated event-times. 
	eventTimes <- timeAtExit[isCase]
	if (any(duplicated(eventTimes))) {  # (or anyDuplicated)
		if (!is.null(precision)) {
			stop("fixTiedEventTimes failed; precision=", precision, 
					" is probably too small")
		} else {	
			nNonUnique <- 
					sum(eventTimes %in% eventTimes[duplicated(eventTimes)])
			stop("there are ", nNonUnique, 
					" non-unique event-times (not allowed), ",
					"but precision=NULL so it was not possible to fix them")
		}
	}
	return(list(
			timeAtEntry=timeAtEntry, timeAtExit=timeAtExit, message=message))
}

################################################################################
# Fix tied values (e.g. event-times) by slightly changing some of them. 

# Details: if there is a pair of tied values, change one of them by 
# epsilonForTies; if there is a triple, increase one of them by epsilonForTies 
# and decrease one; etc.

# NB perturbTiedValues assumes that the values can be compared for equality 
# correctly, using "unique" and "==", even if they are floating-point numbers. 
# If this is not the case, fix them so that they can or change 
# perturbTiedValues to use match(values, unique(values)).  

perturbTiedValues <- function(values, epsilonForTies) {
	numberOfChanges <- 0
	uniqueValues <- unique(values)
	
	for (v in uniqueValues) {
		elements <- which(values == v)
		if (length(elements)==0) 
			stop("INTERNAL ERROR: unique(values) did not work properly")
		numberOfElementsToChange <- length(elements) - 1
		if (numberOfElementsToChange == 0) next
		
		# Example: if numberOfElementsToChange = 3 then amountsToChangeBy = 
		# epsilonForTies * c(1/2, -1/2, 1, -1), or the negative of that (and 
		# because 3 is odd the last element of this will not be used). 
		denom <- ceiling(numberOfElementsToChange / 2)
		posAndNegIncrements <- as.vector(outer(sample(c(-1,1)), 1:denom))
		amountsToChangeBy <- epsilonForTies * posAndNegIncrements / denom 
		
		#cat("Dealing with duplicated value", v, ": elements=", elements, 
		#		" amounts=",amountsToChangeBy[1:numberOfElementsToChange],"\n")
		elements <- sample(elements)  # shuffle (fails if argument length is 1)
		for (i in 1:numberOfElementsToChange) {
			row <- elements[i]
			values[row] <- values[row] + amountsToChangeBy[i]
		}
		numberOfChanges <- numberOfChanges + numberOfElementsToChange
	}
	
	# Make the message and return that and the values. 
	if (numberOfChanges > 0) {
		message <- paste0(numberOfChanges, " of ", length(values), 
				" discretized event-times were changed by up to ", 
				epsilonForTies, " to deal\n with ties.")
	} else { 
		message <- ""
	}
	return(list(values=values, message=message))
}

################################################################################

