################################################################################
# cchs.R

# Functions defined in each R file:
#   cchs.R
#      cchs
#   argProcessingSubroutines.R 
#      checkCchsArgumentsNotMissing
#      checkCchsArguments
#      getExtraVariables
#      getCoxphControl
#      checkSubcohortStrataNotEmpty
#   cchs-class.R
#      print.cchs
#      summary.cchs
#   tiedEventTimeSubroutines.R
#      fixTiedEventTimes
#      perturbTiedValues
#   rowManipulationSubroutines.R
#      doSplitting
#      splitSubcohortCasesJustBeforeEvent
#      dropRowsThatAreNeverAtRisk
#   samplingFractionSubroutines.R
#      checkSamplingFractions
#      getSamplingFractions
#      checkCohortStratumSizes
#      adjustSamplingFractions
#   utils.R
#      allEqual
#      minDifference
#      hasUnusedLevels
#      removeNAsFromVector
#      setS3class

# R style:
# - boolean variables (e.g. isCase, inSubcohort) should be logicals (FALSE/
#   TRUE), not 0/1.  
# - to get a column called b from a data-frame called a, use a$b instead of 
#   a[["b"]], and use $ or [[ in preference to a[,b] 
# - stop("INTERNAL ERROR: ...") means something is wrong with the code, not  
#   with what the user did

# been renamed. I think these are here because I accidentally deleted 
# borganIII folder and then recreated it by copying from C:\MyDocs_backupOnC, 
# and that contains excess files due to having been created using xcopy. 

################################################################################
# Main function to calculate Borgan's Estimator III. 

cchs <- function(formula, data=sys.parent(), inSubcohort, stratum, 
			samplingFractions, cohortStratumSizes, precision=NULL, 
			swap=TRUE, dropNeverAtRiskRows=TRUE, dropSubcohEventsDfbeta=FALSE, 
			adjustSampFracIfAnyNAs=FALSE, keepAllCoxphElements=FALSE, 
			confidenceLevel=0.95, verbose=FALSE, annotateErrors=TRUE, 
			coxphControl, ...) {
	# Note: the "0.95" is hard-coded in cchs-class.R and cchs.Rd as well. 
	
	############################################################################
	# Deal with some of the arguments, etc.  
	originalCall <- match.call() 

	# Check that the required arguments are not missing.
	checkCchsArgumentsNotMissing(formula, inSubcohort, stratum)
	
	# Get the extra variables (inSubcohort, stratum, etc.), which can be columns 
	# of "data", standard variables, or expressions. 
	result <- getExtraVariables(data, originalCall, 
			samplingFractions, cohortStratumSizes, verbose)
	if (length(result) > 0) 
		for (varName in names(result)) 
			assign(varName, result[[varName]])
	
	# Check several of the arguments, and make inSubcohort be boolean/logical. 
	# (Checking class(formula) must be done before identifying NAs.)
	checkCchsArguments(formula, data, inSubcohort, stratum, confidenceLevel) 
	inSubcohort <- as.logical(inSubcohort) 
	
	############################################################################
	# Find all rows that have NAs in the columns of data that appear in formula, 
	# including the event/censoring variable. (This throws an error in certain 
	# situations, so add appropriate extra text to those error-messages, if 
	# annotateErrors=TRUE.) See also extra comment 1 below.
	if (annotateErrors) {
		tryCatch(
			modelFrame <- model.frame(formula, data, na.action=na.omit),
			error=function(e) {
				if (e$call[[1]] == "Surv") {
					errorReason <- "the Surv object being illegal"
				} else if (grepl("^object .* not found$", e$message)) {
					errorReason <- "a variable in the model-formula\nnot existing"
				} else if (grepl("^variable lengths differ", e$message)) {
					errorReason <- 
						"variables in the model-formula having\ndifferent lengths"
				} else {
					errorReason <- "a problem with the model-formula"
				}
				stop(e$message, "\n[This error is caused by ", errorReason, ".]")
			}
		)
	} else {
		modelFrame <- model.frame(formula, data, na.action=na.omit)
	}
	rowsToDrop <- attr(modelFrame, "na.action")
	if (is.null(rowsToDrop)) {
		rowsToDrop <- numeric(0)
	} else { 
		rowsToDrop <- as.vector(rowsToDrop)  # these are all distinct
	}

	# Find rows that have NAs in inSubcohort or stratum and include these too.  
	rowsToDrop <- union(rowsToDrop, which(is.na(inSubcohort)))
	rowsToDrop <- union(rowsToDrop, which(is.na(stratum)))
	
	# Create rowsToKeep, nRowsToKeep, etc. This has to be done before 
	# adjustSamplingFractions. 
	n <- length(inSubcohort)
	rowsToKeep <- setdiff(1:n, rowsToDrop)
	nRowsToDrop <- length(rowsToDrop)
	droppedRowsMessage <- ifelse(nRowsToDrop==0, "", 
			paste0(nRowsToDrop, " observation", ifelse(nRowsToDrop==1, 
			" was", "s were"), " dropped because of NAs."))
	
	# Make stratum a factor. 
	if (!is.factor(stratum)) {
		stratum <- tryCatch(
			as.factor(stratum),
			error=function(e) stop("stratum cannot be converted to a factor")
		)
	}
	
	# Drop any unused levels from stratum. See extra comment 2 below. 
	stratum <- droplevels(stratum) 

	# Check that all strata have at least one row in the subcohort. This will  
	# be assumed by checkCohortStratumSizes. 
	checkSubcohortStrataNotEmpty(stratum, inSubcohort,
		message="before dropping NA-containing rows")

	# Get coxphControl, either from the coxphControl argument or from "...". 
	coxphControl <- getCoxphControl(coxphControl, ...)

	# If verbose, show the first few rows. 
	#if (verbose && is.data.frame(data)) {
	#	cat("## First few rows of data:\n")	
	#	print(head(data))
	#}

	############################################################################
	# Deal with the sampling fractions. 

	# Check samplingFractions or cohortStratumSizes, whichever was given. 
	samplingFractionInfo <- checkSamplingFractions(samplingFractions, 
			cohortStratumSizes, stratum, inSubcohort, verbose)
	
	# Find samplingFractions as a long vector. 
	samplingFractions <- getSamplingFractions(samplingFractions, 
			cohortStratumSizes, stratum, inSubcohort, samplingFractionInfo)
	
	# Check that all cohortStratumSizes >= case-cohort stratum sizes. 
	checkCohortStratumSizes(samplingFractions, cohortStratumSizes, 
			stratum, inSubcohort, verbose=FALSE)  # (or verbose=verbose)
	
	# Adjust samplingFractions, if the user requested it and there were any NAs. 
	if (adjustSampFracIfAnyNAs && nRowsToDrop > 0)
		samplingFractions <- adjustSamplingFractions(samplingFractions, 
				stratum, inSubcohort, rowsToKeep, verbose)  
	
	############################################################################
	# Remake modelFrame, which contains the response and the model matrix, and 
	# drop the NA-containing rows from it. 
	modelFrame <- model.frame(formula, data, na.action=na.pass)
	modelFrame <- modelFrame[rowsToKeep,]
	# The previous line drops the appropriate rows from all attributes of 
	# modelFrame, so the lines with model.matrix and model.response will work.
	if (nrow(modelFrame) == 0) stop("no observations")
	
	# Drop the appropriate rows from other variables. 
	stratum <- stratum[rowsToKeep]
	inSubcohort <- inSubcohort[rowsToKeep]
	samplingFractions <- samplingFractions[rowsToKeep]
	
	# Drop unused levels from stratum again. See extra comment 2 below.
	stratum <- droplevels(stratum) 
	
	# Check that all strata still have at least one row in the subcohort.
	checkSubcohortStrataNotEmpty(stratum, inSubcohort, 
			message="after dropping NA-containing rows")
	
	# Make the model-matrix from the model-frame. 
	modelMatrix <- model.matrix(attr(modelFrame, "terms"), modelFrame)
	modelMatrix <- modelMatrix[, -1, drop=FALSE]  # remove "(Intercept)" column 
	# (the "-1" is probably safer than explicitly using the column's name) 
	
	# Deal with special terms in the model-formula. "cluster" and "strata" are 
	# not allowed but "offset" is. (coxph also has "tt", which is similar but  
	# is passed as an argument to coxph, not in the model-formula.) 
	prohibitedSpecialNames <-  c("cluster", "strata")
	prohibitedSpecialTermsList <- 
			attr(terms(formula, specials=prohibitedSpecialNames), "specials")  
			# "terms(modelFrame, specials=prohibitedSpecialNames)" does not work
	for (i in prohibitedSpecialNames)
		if (!is.null(prohibitedSpecialTermsList[[i]]))
			stop("formula is not allowed to contain these special terms: ", 
					paste(prohibitedSpecialNames, collapse=", "))

	############################################################################
	# Get the Surv object and from that get entry-time, exit-time, and event/
	# censoring variable. (Surv coerces the event/censoring variable to 0/1.) 
	survObject <- model.response(modelFrame)
	if (!inherits(survObject, "Surv")) 
		stop("the response has to be a Surv object")
	censoringType <- attr(survObject, "type")
	if (censoringType != "right" && censoringType != "counting")
		stop("censoring type \"", censoringType, "\" is not allowed")
	if (censoringType == "right") {
		if (!identical(colnames(survObject), c("time", "status")))
			stop("INTERNAL ERROR: unable to process Surv object")
		timeAtEntry <- rep(0, nrow(survObject))
		timeAtExit <- survObject[,"time"]  # survObject$time does not work
		isCase <- survObject[,"status"]
	} else {
		if (!identical(colnames(survObject), c("start", "stop", "status")))
			stop("INTERNAL ERROR: unable to process Surv object")
		timeAtEntry <- survObject[,"start"]
		timeAtExit <- survObject[,"stop"]
		isCase <- survObject[,"status"]
	}
	# Note on Surv types: "counting" is intended for the scenario where each 
	# patient has multiple rows, one for each time-interval. But it is also 
	# suitable for the scenario where each patient has a non-zero entry-time as 
	# well as the exit-time (the exit-time is when they are censored or have 
	# the event). The former scenario is "duration as timescale" and the latter 
	# scenario is "age as timescale". 
	# Names of corresponding variables in cch: 
	# cens=isCase, texit=timeAtExit, tentry=timeAtEntry, X=modelMatrix
	
	# Convert isCase to logical and check there are no non-subcohort non-cases. 
	isCase <- as.logical(isCase) 
	isNonsubcohortNoncase <- !inSubcohort & !isCase
	if (any(isNonsubcohortNoncase)) 
		stop("there are ", sum(isNonsubcohortNoncase), 
				" non-subcohort non-cases")

	# Store numbers of cases etc., which will be put in the final result object. 
	n <- nrow(modelMatrix)
	nEachStatus <- c(
			subcohortNoncases=sum(inSubcohort & !isCase), 
			subcohortCases=sum(inSubcohort & isCase), 
			nonsubcohortCases=sum(!inSubcohort)
	)
	
	############################################################################
	# Deal with tied event-times (and check all times using precision).
	fixedTimes <- fixTiedEventTimes(timeAtEntry, timeAtExit, isCase, precision)
	timeAtEntry <- fixedTimes$timeAtEntry
	timeAtExit <- fixedTimes$timeAtExit
	tiedTimesMessage <- fixedTimes$message 
	
	# Having dealt with tied event-times, now define epsilon, which will be used 
	# (a) in doSplitting, (b) for changing entry-times for non-subcohort cases, 
	# and (c) in splitSubcohortCasesJustBeforeEvent. 
	entryAndExitTimes <- c(timeAtEntry, timeAtExit)
	epsilon <- minDifference(entryAndExitTimes) / 10
	if (verbose) cat("## epsilon=", epsilon, "\n", sep="")
	if (epsilon == 0) 
		stop("epsilon has underflowed to zero, because some entry and \n",
				"exit times were too close together")
	
	# Define and possibly display stratumSizes and subcohortStratumSizes. 
	stratumSizes <- tapply(stratum, stratum, length) 
	subcohortStratumSizes <- sapply(X=levels(stratum), 
			FUN=function(x) sum(inSubcohort[stratum==x])) 
	if (verbose) {
		sf <- round(samplingFractions[match(levels(stratum), stratum)], 6)
		cat("\n## stratumSizes=", paste(stratumSizes,collapse=" "), 
		"\n## subcohortStratumSizes=", paste(subcohortStratumSizes,collapse=" "), 
		"\n## samplingFractions=", paste(sf,collapse=" "),	"\n\n", sep="")
	}
	
	############################################################################
	# Store the column-names of modelMatrix, and replace the column-names by 
	# dummy names so that they are different from all the extra variables' names 
	modelMatrixColumnNames <- colnames(modelMatrix)
	colnames(modelMatrix) <- paste0("dummyName", 1:ncol(modelMatrix))
	
	# Rename modelMatrix as modelMatrixPlus and make it be a data-frame that  
	# also contains the "extra variables" such as inSubcohort and timeAtEntry. 
	# See extra comments 3 and 4 below. 
	modelMatrixPlus <- modelMatrix
	rm(modelMatrix) 
	extraVariables <- c("timeAtEntry", "timeAtExit", "isCase", "inSubcohort", 
			"stratum", "samplingFractions", "id") 
	id <- 1:n  # this will be used in the variance calculations
	modelMatrixPlus <- as.data.frame(modelMatrixPlus)
	for (varName in extraVariables) 
		modelMatrixPlus[[varName]] <- get(varName) 

	# Do the splitting (if swap=TRUE). 
	if (swap) 
		modelMatrixPlus <- doSplitting(modelMatrixPlus, epsilon,verbose=verbose)

	# Change the entry-time for the non-subcohort cases. 
	nonSubcohortRows <- which(!modelMatrixPlus$inSubcohort)
	modelMatrixPlus$timeAtEntry[nonSubcohortRows] <- 
			modelMatrixPlus$timeAtExit[nonSubcohortRows] - epsilon
	
	# Remove any rows that have timeAtEntry=timeAtExit.  
	modelMatrixPlus <- modelMatrixPlus[
			modelMatrixPlus$timeAtEntry != modelMatrixPlus$timeAtExit,] 
	
	# Create modelMatrixPlus$useForDfbeta. This indicates which rows of  
	# modelMatrixPlus will be used in the calculations with the dfbeta  
	# residuals to calculate the variances and CIs. 
	modelMatrixPlus$useForDfbeta <- modelMatrixPlus$inSubcohort
	
	# Make two adjustments if the user requested them. The point at which these 
	# subroutines are called is important. For example,
	# dropRowsThatAreNeverAtRisk assumes that the splitting has been done.
	if (dropSubcohEventsDfbeta)
		modelMatrixPlus <- 
				splitSubcohortCasesJustBeforeEvent(modelMatrixPlus, epsilon)
	if (dropNeverAtRiskRows)
		modelMatrixPlus <- 
				dropRowsThatAreNeverAtRisk(modelMatrixPlus, verbose=FALSE)
		# (verbose=TRUE produces a huge amount of output)
	
	#######################################################################
	# Avoid R CMD check problem (http://stackoverflow.com/q/9439256/1310503):  
	useForDfbeta <- NULL 

	# Remake modelMatrix from modelMatrixPlus, by removing the extra columns, &
	# make those be normal variables again. (After doing this, do not reorder or
	# manipulate the rows of modelMatrix or the elements of inSubcohort etc.) 
	# First, get the extra variables. 
	varNames <- c(extraVariables, "useForDfbeta") 
	for (varName in varNames)  
		assign(varName, modelMatrixPlus[[varName]]) 
	# Remove the extra columns and rename it as modelMatrix. 
	modelMatrix <- modelMatrixPlus[,setdiff(names(modelMatrixPlus), varNames)]
	rm(modelMatrixPlus)  # see extra comment 4 again
	# Convert to matrix (necessary for coxph). 
	modelMatrix <- as.matrix(modelMatrix) 
	# Put the column-names back on modelMatrix. 
	colnames(modelMatrix) <- modelMatrixColumnNames
	
	# Fit the Cox model. coxph sometimes throws errors and warnings, so append 
	# extra text to these to make it clear to the user of cchs where they come 
	# from. 
	# (Intercepting and annotating of errors and warnings, which is done if 
	# annotateErrors=TRUE, works as follows. 
	# If coxph throws an error, append extra text and throw the error.
	# If it throws a warning, append extra text and throw the warning, but still
	# store the result of coxph in result. This code for intercepting errors 
	# and warnings is based on http://stackoverflow.com/a/4952908/1310503.) 
	if (annotateErrors) {
		result <- withCallingHandlers(
			tryCatch(
				coxph(Surv(timeAtEntry, timeAtExit, isCase) ~ modelMatrix + 
					offset(-log(samplingFractions)), control=coxphControl),
				error=function(e) stop(e$message, 
				"\n[This error was thrown by coxph, which was called by cchs.]")
			), 
			warning=function(w) {
				warning(conditionMessage(w), 
				"\n[This warning was thrown by coxph, which was called by cchs.]")
				invokeRestart("muffleWarning")
			}
		)
	} else {
		result <- coxph(Surv(timeAtEntry, timeAtExit, isCase) ~ modelMatrix + 
				offset(-log(samplingFractions)), control=coxphControl)
	}
	
	# Note on the call to coxph: when a matrix/data-frame is used in the RHS of 
	# a formula, this is taken to mean the columns separated by + (see ?formula) 

	# Fix the names on some elements of result. (This needs to be done because 
	# when a matrix is used in the RHS of a model formula, its name gets 
	# prefixed to the coefficient names.)
	names(result$coefficients) <- modelMatrixColumnNames
	dimnames(result$var) <- list(modelMatrixColumnNames, modelMatrixColumnNames)
	
	############################################################################
	# Calculate the variance matrix, using equation (5) in Langholz & Jiao 2007.

	dfbetaUncombined <- residuals(result, type="dfbeta")
	if (is.matrix(dfbetaUncombined)) {
		dfbetaUncombined <- dfbetaUncombined[useForDfbeta, , drop=FALSE]
	} else {  # it must be a vector
		dfbetaUncombined <- matrix(dfbetaUncombined[useForDfbeta], ncol=1)
	} 
	dfbetaCombined <- aggregate(dfbetaUncombined, 
			list(stra=stratum[useForDfbeta], id=id[useForDfbeta]), sum)
	dfbetaStrata <- dfbetaCombined$stra 
	# dfbetaStrata[i] is now the stratum of the ith person in dfbetaCombined
	dfbetaCombined <- dfbetaCombined[, !(names(dfbetaCombined) %in% 
			c("stra","id")), drop=FALSE]  # drop the "stra" & "id" columns 

	nParameters <- ncol(dfbetaUncombined) 
	samplingFractions <- samplingFractions[match(levels(stratum), stratum)] 
	# (both subcohortStratumSizes and samplingFractions have the same order as 
	# levels(stratum))

	stratumWeights <- subcohortStratumSizes * (1 - samplingFractions)
	varAdjustment <- matrix(0, nParameters, nParameters) 
	for (i in levels(stratum)) {
		result$var <- result$var + stratumWeights[i] * 
				var(dfbetaCombined[dfbetaStrata==i,])
	}
	
	############################################################################
	# Create result$coeffsTable, which contains the hazard ratios, confidence 
	# intervals, p-values, etc. This is largely based on code from 
	# survival:::summary.coxph. See extra comment 5 below. 
	beta <- result$coefficients
	se <- sqrt(diag(result$var))
	z <- qnorm((1 + confidenceLevel) / 2)  # e.g. z=1.96
	result$coeffsTable <- cbind(exp(beta), exp(beta-z*se), exp(beta+z*se), 
			1-pchisq((beta/se)^2, 1), beta, se)		
	colnames(result$coeffsTable) <- 
			c("HR", "CIlower", "CIupper", "p", "logHR", "SElogHR")
	# The p-value is for the Wald test. See Collett 1994/2003, "Modelling 
	# Survival Data in Medical Research", second para of section 3.4, p67/69. 
	
	# Fix other contents of result. 
	result$call <- originalCall
	result$n <- n 
	result$nEachStatus <- nEachStatus
	result$nevent <- n - result$nEachStatus["subcohortNoncases"] 
	result$nStrata <- nlevels(stratum)  # (NB empty strata have been dropped)
	
	# Make result$message, which is about rows that are dropped because of NAs 
	# and event-times that have been changed to avoid ties. 
	message <- droppedRowsMessage 
	if (message != "" && tiedTimesMessage != "") 
		 message <- paste0(message, "\n")
	if (tiedTimesMessage != "")
		message <- paste0(message, tiedTimesMessage)
	result$message <- message
		
	# Store confidenceLevel, if it is not 0.95. 
	if (confidenceLevel != 0.95) 
		result$confidenceLevel <- confidenceLevel
	
	# If !keepAllCoxphElements, then drop the appropriate elements from result. 
	varsToKeep <- c("coefficients","var","iter","n","nevent","nStrata","call",
			"coeffsTable","nEachStatus","message","confidenceLevel")
	if (!keepAllCoxphElements) 
		for (varName in names(result))
			if (!(varName %in% varsToKeep))
				result[varName] <- NULL
	
	# Set the S3 class of result to "cchs" and return it.
	result <- setS3class(result, "cchs")  # works in S as well as R
	return(result)
}

################################################################################
# Extra comments 

# 1. modelFrame <- model.frame(...) is done twice, the first time to detect 
# "rows" that would be used in the model but contain NAs, the second time to  
# create the model-frame that will be manipulated and passed to coxph. The 
# reason for not doing it once to do both these things is as follows. The best 
# way to find what formula-variables have NAs is to run model.frame with 
# na.action=na.omit, and doing this also removes those rows from the resulting 
# model-frame. (If you do model.frame with na.pass, then attr(modelFrame, 
# "na.action") does not get created.) When making the model-frame that will be  
# manipulated and passed to coxph, it is cleaner and safer to use "na.pass" 
# (which means don't drop the NA-containing rows) and then remove all the 
# NA-containing rows from modelMatrix, inSubcohort, etc. in one fell swoop. 
# The NA-containing rows are identified by rowsWithNA. 
# 
# There is an alternative way of dealing with NAs and which rows should be used.
# Find which rows are dropped by coxph by including cluster(1:nrow(data)) in 
# the formula and using model=TRUE, and then look at the "cluster(1:nrow(data))"
# column of the coxph object. (It is then also necessary to use result$naive.var
# instead of result$var when calculating the variance matrix.) 
# However, creating modelFrame twice has various advantages: it is easier to 
# catch errors in the Surv object or model-formula and easier to identify 
# which elements of inSubcohort etc. need to be kept. 

# 2. The assumption that all the levels of stratum appear in the data, i.e. that 
# stratum has no unused levels like you get from as.factor(4:6)[-2], is made by:
# - the functions in samplingFractionSubroutines.R 
# - the loop that finally calculates the correct value of result$var
# - for (stra in levels(stratum)) ...
# and maybe other code after the first droplevels as well. Making this 
# assumption makes those other pieces of code much more readable. For it to 
# hold, it is necessary to do droplevels twice, firstly before 
# checkSamplingFractions and secondly after the rowsToDrop are removed. 

# 3. Reasons for making modelMatrixPlus a data-frame rather than a matrix:
# - data-frames can store different types of variables, so stratum (a 
#   factor), inSubcohort (a boolean/logical), etc. will retain their types 
#   when they are put into modelMatrixPlus and later removed from it;
# - it allows you to write m$b in doSplitting and other subroutines, which is
#   more concise than m[,b], and names(m) rather than colnames(m);
# - having all these things in a single data-frame makes it easier to do 
#   the splitting and other manipulations. 

# 4. I think "a <- b; rm(b)" is the best way to have meaningful variable-names 
# without wasting a lot of memory. See also 
# http://stackoverflow.com/a/2717853/1310503, and search for "lazy evaluation" 
# in R Language Definition.

# 5. Notes on the internal code in survival:::summary.coxph: se is calculated 
# from myModel$var and the CIs are calculated using se. In rval$coefficients, 
# se is called "robust se", but in cchs this would be a misnomer because cchs 
# stores the asymptotic variance in myModel$var. 

################################################################################

