################################################################################
# Subroutines to do manipulation of the augmented model-matrix 
# (modelMatrixPlus), before the calculation of the estimator. 

################################################################################
# Subroutine to do the splitting. 

doSplitting <- function(modelMatrixPlus, epsilon, verbose=FALSE) { 
	# Avoid R CMD check problem (http://stackoverflow.com/q/9439256/1310503):  
	stratum <- inSubcohort <- timeAtExit <- timeAtEntry <- NULL
	
	# Do some checking. 
	mmp <- modelMatrixPlus
	if (!is.data.frame(mmp))
		stop("the modelMatrixPlus that is passed to doSplitting has to be a ",
				"data-frame")
	if (any(duplicated(mmp$id)))
		stop("the modelMatrixPlus that is passed to doSplitting has to be ",
				"an\naugmented model matrix in which all rows correspond ",
				"to different\nrows of the original data") 
	# Check that mmp contains all the required extra variables. 
	extraVariables <- c("stratum", "isCase", "inSubcohort", "timeAtEntry", 
			"timeAtExit")
	if (!all(extraVariables %in% colnames(mmp)))
		stop("modelMatrixPlus has to have the following columns:\n", 
				paste(extraVariables, collapse=" "))
	
	# Make it possible to use "inSubcohort" etc., because this is more concise. 
	# Do this for all columns whose names don't start with "dummyName".  
	for (varName in grep("^dummyName", names(mmp), value=TRUE, invert=TRUE)) 
		assign(varName, mmp[[varName]])
	if (!is.factor(stratum)) 
		stop("modelMatrixPlus$stratum has to be a factor")
	nStrata <- nlevels(stratum)
	# Other possible checks: timeAtEntry < timeAtExit, etc.
	
	# Decide which row in each stratum should be the one that gets split. 
	# (Record this using row number. This assumes you are not going to re-order 
	# the rows between choosing which rows to split and doing the splitting.)
	rowToSplit <- numeric(nStrata) 
	names(rowToSplit) <- levels(stratum)
	# rowToSplit["x"] will be the row to split for non-subcohort cases in 
	# stratum x
	if (verbose) 
		cat("######### doSplitting #########\nChoose rows to split:\n")
	for (i in levels(stratum)) {
		possibleRows <- which(inSubcohort & stratum==i)  
		# possibleRows is the rows who are in the subcohort and the same stratum
		if (length(possibleRows) == 0) {  # (this happens if there is no one in 
			rowToSplit[i] <- NA            # the subcohort and the same stratum)
		} else if (length(possibleRows) == 1) {  # (see "undesired behaviour" in
			rowToSplit[i] <- possibleRows        # ?sample)
		} else { 
			rowToSplit[i] <- sample(possibleRows, size=1)
		}
		if (verbose) 
			cat("  Row to split in stratum ", i, ": row ", rowToSplit[i], 
				" (possibleRows=", paste(possibleRows,collapse=","), ")\n", 
				sep="")
	}
	if (verbose) cat("\n")
	
	# Decide where each rowToSplit needs to be split. 
	if (verbose) 
		cat("For each non-subcohort case, work out what split to do:\n")
	splitTimes <- vector("list", nStrata)  # make empty list
	names(splitTimes) <- levels(stratum)
	# splitTimes$str3 will be the times at which row number rowToSplit$str3 
	# needs to be split
	for (i in which(!inSubcohort)) { 
		# Deal with the non-subcohort case in row i. 
		thisStratum <- stratum[i]
		if (is.na(rowToSplit[thisStratum])) next 
		eventTime <- timeAtExit[i]
		if (verbose) cat("  Non-subcohort case in row ", i, " (stratum ", 
				as.character(thisStratum), "): rowToSplit=", 
				rowToSplit[thisStratum], " will be split at eventTime=", 
				eventTime, " (if it is at risk then)\n", sep="")
		# If this stratum's rowToSplit is at risk at eventTime, record that a 
		# split needs to be made.
		if (timeAtEntry[rowToSplit[thisStratum]] < eventTime && 
				eventTime <= timeAtExit[rowToSplit[thisStratum]])
			splitTimes[[thisStratum]] <- c(splitTimes[[thisStratum]], eventTime)
	}
	
	# Do the splitting. 
	if (verbose) cat("\nDo splits:\n")
	for (i in levels(stratum)) { 
		if (is.null(splitTimes[[i]])) next  # there are no times to split at
		sortedSplitTimes <- sort(splitTimes[[i]])
		originalRow <- rowToSplit[i] 
		if (verbose) 
			cat("  Stratum ", i, ": splitting row ", originalRow, " at times ", 
					paste(sortedSplitTimes,collapse=" "), "\n", sep="")
		for (j in 1:length(sortedSplitTimes)) {
			newRow <- nrow(mmp) + 1
			mmp[newRow,] <- mmp[originalRow,]  # create the new row
			mmp$timeAtEntry[originalRow] <- sortedSplitTimes[j] 
			mmp$timeAtExit[newRow] <- sortedSplitTimes[j] - epsilon  
			# The original row now contains the time from split-time onwards,  
			# and the new row contains the time before the split-time. 
			mmp$isCase[newRow] <- FALSE  # the new row does not have an event
		}
	}
	if (verbose) cat("###############################\n\n")
	
	return(mmp)
}

################################################################################
# Subroutines for two adjustments to do with dfbeta residuals. See 
# dropSubcohEventsDfbeta and dropNeverAtRiskRows in the Arguments section of 
# ?cchs. 

# 1. For subcohort cases, split off the time at the very end, so that these 
# short time-segments can be excluded when calculating dfbeta. This needs to be 
# done after the splitting.
splitSubcohortCasesJustBeforeEvent <- function(modelMatrixPlus, epsilon) {
	# (No checking, since this is a fairly short function.)
	mmp <- modelMatrixPlus
	for (originalRow in which(mmp$inSubcohort & mmp$isCase)) { 
		mmp <- rbind(mmp, mmp[originalRow,])  # add a copy of originalRow
		newRow <- nrow(mmp)
		mmp$timeAtExit[originalRow] <- 
				mmp$timeAtExit[originalRow] - epsilon 
		mmp$timeAtEntry[newRow] <- mmp$timeAtExit[newRow] - epsilon
		mmp$isCase[originalRow] <- FALSE  # the original row is not a case
		mmp$useForDfbeta[newRow] <- FALSE
	}
	return(mmp)
}

# 2. Drop rows that are never at risk, so that they will not be used by coxph 
# (and obviously they will not be used in the dfbeta calculations either). 
dropRowsThatAreNeverAtRisk <- function(modelMatrixPlus, verbose=FALSE) {
	# (No checking, since this is a short function.)
	keepTheseRows <- logical(nrow(modelMatrixPlus))
	eventTimes <- modelMatrixPlus$timeAtExit[which(modelMatrixPlus$isCase)]
	if (verbose) 
		cat("dropRowsThatAreNeverAtRisk  eventTimes:", eventTimes, "\n")
	for (i in 1:nrow(modelMatrixPlus)) {
		# Keep this row if its at-risk time contains any of the event times:
		keepTheseRows[i] <- any(modelMatrixPlus$timeAtEntry[i] < eventTimes & 
				eventTimes <= modelMatrixPlus$timeAtExit[i])
		if (verbose) 
			cat(" row", i, "is at risk at eventTimes?", 
					(modelMatrixPlus$timeAtEntry[i] < eventTimes & 
					eventTimes <= modelMatrixPlus$timeAtExit[i]), "\n")
	}
	return(modelMatrixPlus[keepTheseRows,])
}

################################################################################
