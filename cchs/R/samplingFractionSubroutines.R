################################################################################
# Subroutines for checking and processing the sampling fractions. 

################################################################################
# Check samplingFractions/cohortStratumSizes in various ways. 

checkSamplingFractions <- function(samplingFractions, cohortStratumSizes, 
		stratum, inSubcohort, verbose=FALSE) {
	# When this function is run, none of the arguments will yet have had the 
	# NA-containing rows (rowsToDrop) removed. 
	
	# Check stratum and inSubcohort are the same length, and define some things. 
	if (length(stratum) != length(inSubcohort))
		stop("INTERNAL ERROR: stratum and inSubcohort have to be the same ",
				"length")
	if (!is.factor(stratum))
		stop("INTERNAL ERROR: stratum has to be a factor")
	if (hasUnusedLevels(stratum))
		stop("INTERNAL ERROR: stratum should not have any unused levels")
	nObservations <- length(stratum)  # = length(inSubcohort)
	nStrata <- nlevels(stratum)
	sfGiven <- !missing(samplingFractions)
	cssGiven <- !missing(cohortStratumSizes)
	
	# Check that exactly one of the two variables was given, and call it arg.
	if (!xor(sfGiven, cssGiven)) 
		stop("exactly one of samplingFractions and cohortStratumSizes must ", 
				"be supplied") 
	argName <- ifelse(sfGiven, "samplingFractions", "cohortStratumSizes")
	arg <- get(argName)  # "arg <- ifelse(sfGiven, ..." does not work	
	
	# Check that arg is legal. 
	# If arg has no names, then it has to have length nObservations, and arg  
	# is "long"; and arg is not allowed to have NAs. 
	# If arg has names, then these have to include the levels of stratum, and 
	# arg is "short"; and arg is only allowed to have NAs in unused elements,  
	# meaning ones whose names do not appear in levels(stratum). 
	
	# First, make an error message to use when arg has an illegal form. 
	if (is.null(names(arg))) {
		namesArgText <- "NULL"
	} else {
		namesArgText <- paste(paste(names(arg)[1:min(6,length(arg))], 
						collapse=" "),	if(length(arg)>6) "...")
	}
	levelsStratumText <- paste(paste(levels(stratum)[1:min(6,nStrata)], 
					collapse=" "), if(nStrata>6) "...")
	illegalArgumentFormMessage <- paste0(argName, " must either have names ",
			"that correspond to the levels of stratum\nand are all distinct, ",
			"or have no names and length equal to the number of ",
			"observations \nnames(", argName, ")=", namesArgText, 
			"\nlevels(stratum)=", levelsStratumText)
	
	# Now check arg in the two cases where it has no names and it has names. 
	# Check the length, names (if any), and NAs. 
	if (is.null(names(arg))) {
		# arg has no names (NAs are not allowed)
		if (length(arg) != nObservations) 
			stop(illegalArgumentFormMessage)
		if (any(is.na(arg)))
			stop(argName, " has no names, so it is not allowed to contain NAs")
		argIsShort <- FALSE
	} else {
		# arg has names (NAs in arg are only allowed in unused elements)
		if (!all(levels(stratum) %in% names(arg)) || 
				any(duplicated(names(arg)))) 
			stop(illegalArgumentFormMessage)
		illegalNaNames <- intersect(names(arg)[is.na(arg)], levels(stratum))
		if (length(illegalNaNames) > 0) 
			stop("elements of ", argName, " whose names match values of ",
					"stratum must not be NA\nIllegal elements: ", 
					paste(illegalNaNames, collapse=" "))
		if (length(arg) > nStrata) 
			cat("These elements of", argName, "will not be used:", 
					setdiff(names(arg), levels(stratum)), "\n")
		arg <- arg[levels(stratum)]  # drop unneeded values (not essential)
		argIsShort <- TRUE
	}
	
	# If arg is long, check it is constant within each stratum. (This checks all 
	# of arg, not just arg[rowsToKeep], because if arg is not consistent within 
	# each stratum then there must be something very wrong with the data.) 
	if (!argIsShort)
		for (stra in levels(stratum))
			if (!allEqual(arg[stratum==stra]))
				stop(argName, " is not constant within stratum ", stra)
	
	# Check that the numerical values of arg are legal. (Check all of arg, not
	# just arg[rowsToKeep].)
	if (sfGiven) {
		if (any(arg <= 0 | arg > 1, na.rm=TRUE))
			stop(argName, " must all be greater than 0 and less than or ",
					"equal to 1")
	} else {
		if (any(arg %% 1 != 0 | arg <= 0, na.rm=TRUE))
			stop(argName, " must all be positive integers")
	}
	
	# Return the two booleans that are needed by getSamplingFractions, which  
	# will be called next. 
	return(list(sfGiven=sfGiven, argIsShort=argIsShort))
	
	# Not done: check that each cohort stratum size times the corresponding 
	# sampling fraction equals an integer.
}

################################################################################
# Get samplingFractions as a "long" vector of length nrow(data), so that it  
# can ultimately be used in the offset term in the call to coxph.

getSamplingFractions <- function(samplingFractions, cohortStratumSizes, 
		stratum, inSubcohort, samplingFractionInfo, verbose=FALSE) {
	
	# Check that stratum has no unused levels. 
	if (hasUnusedLevels(stratum))
		stop("INTERNAL ERROR: stratum should not have any unused levels")
	
	# stratum is allowed to contain NAs in this function, which is why  
	# na.rm=TRUE is needed below. 
	
	# samplingFractionInfo is a list of two booleans. (Don't bother checking!)
	sfGiven <- samplingFractionInfo$sfGiven
	argIsShort <- samplingFractionInfo$argIsShort
	
	# There are four possible situations, corresponding to the two booleans:
	
	if (sfGiven) {  # samplingFractions was given
		if (argIsShort) {
			return(samplingFractions[match(stratum, names(samplingFractions))])
		} else {
			return(samplingFractions)
		}
		
	} else {  # cohortStratumSizes was given
		# Calculate and check subcohortStratumSizes:
		subcohortStratumSizes <- sapply(X=levels(stratum), 
				FUN=function(x) sum(inSubcohort[stratum==x], na.rm=TRUE))
		if (sum(subcohortStratumSizes) == 0) 
			stop("subcohortStratumSizes are all zero")
		
		# Make cohortStratumSizes be short and have the right order:
		if (argIsShort) {  # make it have the right order
			cohortStratumSizes <- cohortStratumSizes[
					match(levels(stratum), names(cohortStratumSizes))]
			# Check that the match in the previous line worked correctly:
			if (!identical(names(subcohortStratumSizes), 
					names(cohortStratumSizes)))
				stop("INTERNAL ERROR: code in getSamplingFractions must ",
						"be wrong") 
		} else {  # make it short  
			cohortStratumSizes <- cohortStratumSizes[match(levels(stratum), 
					stratum)]
		}
		
		# Calculate the sampling fraction for each stratum:
		samplingFractions <- subcohortStratumSizes / cohortStratumSizes
		if (verbose) {
			cat("## getSamplingFractions calculations:\n")
			for (varName in c("subcohortStratumSizes", "cohortStratumSizes", 
					"samplingFractions")) {
				cat("\n  ", varName, ":\n", sep="")
				print(get(varName))
			}
			cat("\n")
		}
		return( samplingFractions[match(stratum, names(samplingFractions))] )
	}
}

################################################################################
# Check that the cohortStratumSizes are greater than the subcohortStratumSizes.  

checkCohortStratumSizes <- function(samplingFractions, cohortStratumSizes, 
		stratum, inSubcohort, verbose=FALSE) {
	# Partly check the arguments that were passed to this function (maybe this  
	# checking is over-the-top). 
	
	if (verbose && !missing(cohortStratumSizes)) 
		cat("## checkCohortStratumSizes:\n   stratum:", stratum, 
				"\n   cohortStratumSizes:", cohortStratumSizes, "\n")
	
	# Check samplingFractions was supplied and has the same length as stratum. 
	if (missing(samplingFractions) || 
			length(samplingFractions) != length(stratum)) 
		stop("INTERNAL ERROR: checkCohortStratumSizes is intended to be ",
				"called after getSamplingFractions, so that ",
				"samplingFractions exists and is \"long\"")
	
	# Check stratum and inSubcohort have the same length. 
	if (length(stratum) != length(inSubcohort))
		stop("stratum and inSubcohort have to be the same length")
	
	# Check samplingFractions and cohortStratumSizes are constant within strata.
	varNames <- "samplingFractions"
	if (!missing(cohortStratumSizes)) {
		varNames <- c(varNames, "cohortStratumSizes")
		if (!is.null(names(cohortStratumSizes))) # it's "short"; make it "long"
			#if (length(cohortStratumSizes) == nlevels(stratum)) #(alternative)
			cohortStratumSizes <- 
					cohortStratumSizes[match(stratum, levels(stratum))]
	}
	for (varName in varNames)
		for (stra in levels(stratum))
			if (!allEqual(removeNAsFromVector(get(varName)[stratum==stra])))
				# removeNAsFromVector is needed in case any stratum=NA
				stop(varName, " is not constant within stratum ", stra)
	
	# Check that stratum has no unused levels. 
	if (hasUnusedLevels(stratum))
		stop("INTERNAL ERROR: stratum should not have any unused levels")
	
	# Not done: if cohortStratumSizes is given, check sf is correct. 
	
	# Calculate subcohortStratumSizes and check it. 
	subcohortStratumSizes <- sapply(X=levels(stratum), 
			FUN=function(x) sum(inSubcohort[stratum==x], na.rm=TRUE))
	if (any(subcohortStratumSizes == 0))
		stop("INTERNAL ERROR: subcohortStratumSizes should all be non-zero")
	
	# If cohortStratumSizes was not supplied, calculate it. 
	if (missing(cohortStratumSizes)) {
		# Make samplingFractions "short". 
		samplingFractions <- samplingFractions[match(levels(stratum), stratum)]
		names(samplingFractions) <- levels(stratum)
		
		# Check names of vectors. 
		if (!identical(names(samplingFractions), names(subcohortStratumSizes)))
			stop("INTERNAL ERROR: names(samplingFractions) should be the ",
					"same as names(subcohortStratumSizes)")
		
		# Calculate cohortStratumSizes. 
		cohortStratumSizes <- round(subcohortStratumSizes / samplingFractions)
		if (verbose) { 
			cat("## checkCohortStratumSizes:\n")
			cat("## subcohortStratumSizes=\n"); print(subcohortStratumSizes) 
			cat("## cohortStratumSizes=\n"); print(cohortStratumSizes) 
		}
		
	} else {
		# cohortStratumSizes was supplied; make it be "short".
		if (length(cohortStratumSizes) == length(stratum)) {
			cohortStratumSizes <- 
					cohortStratumSizes[match(levels(stratum),stratum)]
			names(cohortStratumSizes) <- levels(stratum)
		}
	}
	
	# Check that cohortStratumSizes is not too small. 
	for (stra in levels(stratum)) { 
		cohortStratumSize <- cohortStratumSizes[stra]
		rowsInStratum <- sum(stratum==stra, na.rm=TRUE)
		if (cohortStratumSize < rowsInStratum) 
			stop("cohortStratumSizes[\"", stra, "\"]=", cohortStratumSize, 
					" but stratum ", stra, " has size ", rowsInStratum, 
					" in the case-cohort data")
	}
}

################################################################################
# Adjust samplingFractions to account for the rows that are being dropped 
# because they had NAs (these are known as the NA-containing rows). 

# (Possibly it is wasteful to reconstruct the short version of 
# cohortStratumSizes here, when the user might have given that to cchs 
# originally. There are other possibly wasteful things too. But it is cleaner 
# to have all the code for adjusting samplingFractions in one place, especially 
# as it might be rare for this function to even get used.) 

adjustSamplingFractions <- function(samplingFractions, stratum, inSubcohort, 
		rowsToKeep, verbose=FALSE) {
	# Check the arguments. 
	if (!allEqual(
			sapply(list(samplingFractions, stratum, inSubcohort), length)))
		stop("INTERNAL ERROR: samplingFractions, stratum, and rowsToKeep ",
				"have to all be the same length")
	if (max(rowsToKeep) > length(samplingFractions))
		stop("INTERNAL ERROR: rowsToKeep is inconsistent with other vectors")
	
	# In this function it is permissible for stratum to have unused levels and 
	# stratum to be NA. 
	
	# For each stratum, and for before and after the NA-containing rows are  
	# removed, find subcohortStratumSizes as a "short" vector. 
	subcohortStratumSizesBefore <- sapply(X=levels(stratum), 
		FUN=function(x) sum(inSubcohort[stratum==x], na.rm=TRUE))
	subcohortStratumSizesAfter <- sapply(X=levels(stratum), 
		FUN=function(x) sum(inSubcohort[rowsToKeep][stratum[rowsToKeep]==x]))
	# (The "sum()" in the previous line looks clumsy, but this is unavoidable.)
	if (any(is.na(subcohortStratumSizesAfter)))
		stop("INTERNAL ERROR: rowsToKeep contains rows with stratum=NA")
	
	# Similarly find cohortStratumSizesBefore and cohortStratumSizesAfter.  
	samplingFractionsBefore <- 
			samplingFractions[match(levels(stratum), stratum)]
	cohortStratumSizesBefore <- 
			subcohortStratumSizesBefore / samplingFractionsBefore
	numberToDropInEachStratum <- c(table(stratum) - table(stratum[rowsToKeep]))
	cohortStratumSizesAfter <- 
			cohortStratumSizesBefore - numberToDropInEachStratum
	
	# Calculate samplingFractionsAfter, and return it as a "long" vector. 
	samplingFractionsAfter <- 
			subcohortStratumSizesAfter / cohortStratumSizesAfter
	
	# Display things if required. 
	if (verbose) {
		cat("\n## adjustSamplingFractions:\n")
		for (varName in c(
				"subcohortStratumSizesBefore", "cohortStratumSizesBefore", 
				"subcohortStratumSizesAfter", "cohortStratumSizesAfter", 
				"samplingFractionsBefore", "samplingFractionsAfter")) { 
			cat("## ", varName, ":\n", sep="")
			print(get(varName)) 
		}
		cat("\n")
	}
	
	# Return samplingFractionsAfter as a "long" vector with no names. 
	result <- samplingFractionsAfter[
			match(stratum, names(samplingFractionsAfter))]
	names(result) <- NULL
	return(result)
}

################################################################################

