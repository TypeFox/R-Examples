################################################################################
# Subroutines for checking and processing the arguments that are passed to cchs. 

################################################################################
# Basic checking of the arguments. 

checkCchsArgumentsNotMissing <- function(formula, inSubcohort, stratum) {
	# These arguments are all required, so check they are not missing. 
	if (missing(formula)) 
		stop("formula has to be supplied")
	if (missing(inSubcohort)) 
		stop("inSubcohort has to be supplied")
	if (missing(stratum)) 
		stop("stratum has to be supplied")
	
	# missing does not work properly with loops. Code that shows this:
	#f <- function(x) {
	#	for (varName in "x")
	#		print(missing(varName))
	#}
	#f()
}

checkCchsArguments <- function(formula, data, inSubcohort, stratum, 
		confidenceLevel) {
	# Check formula is a formula. 
	if (class(formula) != "formula")
		stop("formula has to be an object of class formula")
	
	# Check data is not an empty data-frame. 
	if (is.data.frame(data) && nrow(data) == 0) 
		stop("data is empty")
	
	# Check the lengths of inSubcohort and stratum. 
	if (is.data.frame(data)) {
		for (varName in c("inSubcohort","stratum"))
			if (length(get(varName)) != nrow(data))
				stop(varName, " has to have length nrow(data)")
	} else {
		if (length(inSubcohort) != length(stratum))
			stop("inSubcohort and stratum have to be the same length")
	}
		
	# Check inSubcohort. 
	if (!all(inSubcohort %in% c(0,1,NA)))  # works if it is boolean/logical too
		stop("inSubcohort values must all be 0, 1, TRUE, or FALSE (or NA)")
	
	# Check stratum. 
	if (length(unique(stratum)) == 1)
		stop("only one stratum; use cch with method=\"Prentice\" instead")
	
	# Check confidenceLevel. 
	if (!is.numeric(confidenceLevel) || length(confidenceLevel) != 1 || 
			confidenceLevel < 0 || confidenceLevel > 1)
		stop("confidenceLevel must be a single numeric in the range [0,1]")
}

################################################################################
# Get inSubcohort and stratum from data, if data is a data-frame and the user 
# specified a variable that is a column of it. Also do this for the one of 
# samplingFractions and cohortStratumSizes that was specified by the user. 
# (It might be possible to do this by making a formula that contains these 
# variables and using model.frame, but the user can provide inSubcohort as
# an explicit expression, for example, so it might not.) 

getExtraVariables <- function(data, originalCall, samplingFractions, 
		cohortStratumSizes, verbose=FALSE) {
	if (verbose) cat("## Dealing with inSubcohort, stratum, etc.:\n")
	
	# Create the set of variables that will be looked for in data. 
	varNames <- c("inSubcohort","stratum")
	if (!missing(samplingFractions)) 
		varNames <- c(varNames, "samplingFractions")
	if (!missing(cohortStratumSizes))
		varNames <- c(varNames, "cohortStratumSizes")
	
	# If data is a data-frame, look for the variables in it. 
	result <- list()
	if (is.data.frame(data)) {
		
		# Example of what the following loop does, with varName="stratum": if 
		# the user typed stratum=xyz, then set stratum to be data$xyz, if 
		# that exists.
		for (varName in varNames) {
			value <- originalCall[[varName]]  # in the example, value=xyz
			valueAsString <- deparse(value)   # in the ex., valueAsString="xyz"
			if (verbose) cat("   ", varName, ": ", sep="")
			if (is.name(value)) {
				# value is a "name" (the name of a variable, not an expression)
				if (valueAsString %in% names(data)) {
					if (verbose) 
						cat("getting", valueAsString, "from the data-frame\n")
					result[[varName]] <- data[[valueAsString]] 
					# in the ex., previous line does result$stratum <- data$xyz
				} else {
					if (verbose) 
						cat(valueAsString," is not in the data-frame\n",sep="")
				}
			} else {
				# value is not a "name"; it might be something like c(0.1, 0.3)
				if (verbose) cat(valueAsString, " is not a \"name\"\n", sep="")
			}
		}
		
	}
	return(result)
}
	
################################################################################
# Get coxphControl, either from the coxphControl argument or from "..." (the 
# latter, if both of these are empty). 

# This also checks coxphControl and converts all warnings to errors. The 
# purpose of checking coxphControl now is to avoid it causing problems later  
# when it is used by coxph. There are many possible problems with  
# coxph.control, some of which it throws as warnings or errors, and some of  
# which are only reported when the object is evaluated. So these are just  
# cursory checks, not thorough. 

getCoxphControl <- function(coxphControl, ...) {
	# Avoid R CMD check problem (http://stackoverflow.com/q/9439256/1310503):  
	coxph.control <- survival::coxph.control
	# Alternative ways that fail: coxph.control <- NULL causes errors below;
	# utils::globalVariables(coxph.control) fails R CMD check
	
	# Check that any extra arguments are legal arguments for coxphControl. 
	extraArgs <- list(...)
	if (length(extraArgs) > 0) {
		coxphArgNames <- names(formals(coxph.control))
		illegalArgIndicators <- is.na(pmatch(names(extraArgs), coxphArgNames))
		if (any(illegalArgIndicators)) {
			illegalArgNames <- names(extraArgs)[illegalArgIndicators]
			stop("illegal argument", if (length(illegalArgNames)>1) "s", 
					" passed to cchs: ", paste(illegalArgNames, collapse=" "))
		}
	}
	
	if (missing(coxphControl)) {
		# Create coxphControl from the extra arguments. 
		tryCatch(
			coxphControl <- coxph.control(...), 
			error=function(e) stop(e$message, 
			"\n[This error was thrown by coxph.control and is due to one ",
			"or more illegal\nvalues of extra arguments (...), which cchs ",
			"passes to coxph.control.]"),
			warning=function(e) stop(e$message, 
			"\n[This error is caused by a warning from coxph.control and ",
			"is due to one or more\nillegal values of extra arguments (...), ",
			"which cchs passes to coxph.control.]")
		)
		
	} else {
		# coxphControl was given, so check that extra arguments were not given. 
		if (length(extraArgs) > 0)
			stop("cchs does not allow both coxphControl and extra arguments ",
					"(...) to\nbe specified")
		# Check coxphControl by evaluating it. 
		tryCatch(
			eval(coxphControl),
			error=function(e) stop(e$message, 
			"\n[This error was thrown when evaluating coxphControl and is ",
			"probably\ndue to use of coxph.control with illegal arguments.]"), 
			warning=function(e) stop(e$message, "\n[This error is caused by ",
			"a warning when evaluating coxphControl and\nis probably due ",
			"to use of coxph.control with illegal arguments.]")
		)
	}
	
	return(coxphControl)
}

################################################################################
# Check that all strata have at least one row in the subcohort.

checkSubcohortStrataNotEmpty <- function(stratum, inSubcohort, message=NULL) {
	if (!is.factor(stratum) || hasUnusedLevels(stratum))
		stop("INTERNAL ERROR: stratum should be a factor with no unused levels")
	if (!is.vector(inSubcohort) || !is.logical(inSubcohort))
		stop("INTERNAL ERROR: inSubcohort should be a vector of logicals")
	if (!is.null(message) && !is.character(message))
		stop("INTERNAL ERROR: message should be either NULL or a string")
	for (stra in levels(stratum)) {
		if (!any(inSubcohort[stratum==stra], na.rm=TRUE)) {
			errorMessage <- paste("subcohort for stratum", stra, "is empty")
			if (!is.null(message)) 
				errorMessage <- paste0(errorMessage, " (", message, ")")
			stop(errorMessage)
		}
	}
}

################################################################################

