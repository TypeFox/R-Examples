################################################################################
# Functions related to the cchs class or classes in general. 

################################################################################
# S3 generic function that will be called when you type cchs(...) or 
# print(cchs(...)). To pass R CMD check, print.cchs must have exactly the same 
# arguments as print.  
print.cchs <- function(x, ...) {
	cat("Call: ")
	print(x$call)
	cat("Number of observations/rows: ", x$n, ", of which\n", sep="")
	statusNumbers <- x$nEachStatus
	cat("        subcohort non-cases:", statusNumbers["subcohortNoncases"], 
			"\n")
	cat("        subcohort cases:    ", statusNumbers["subcohortCases"], 
			"\n")
	cat("        non-subcohort cases:", statusNumbers["nonsubcohortCases"], 
			"\n")
	if (x$message != "") cat(x$message, "\n") 
	cat("Number of strata for subcohort-selection:", x$nStrata, "\n")
	confidenceLevel <- ifelse("confidenceLevel" %in% names(x), 
			x$confidenceLevel, 0.95)
	cat("Coefficients (HR = hazard ratio, CI = ", confidenceLevel * 100, 
			"% confidence interval, SE = standard error):\n", sep="")
	print(x$coeffsTable)
	invisible(x)
}

# S3 generic function that will be called when you type summary(cchs(...)).  
# The purpose of this is just to prevent summary((cchs(...)) from giving 
# useless output. It returns the object unchanged and tells the user what to 
# type to get the information they probably want. 
summary.cchs <- function(object, ...) {
	cat("summary(cchsObject) does no extra calculations and just returns",
			"cchsObject.\nTo get a summary of cchsObject, type cchsObject.\n") 
	#setS3class(object, "summary.cchs")
	return(object)
}

################################################################################
