# CheckMarginals #############
##################
# Returns an ovariable with a marginal devised from the data and upstream variable marginals. 
# Marginal values for data could be stored into the database.
# Assumes that all depended upon variables are in memory, as should be the case.

CheckMarginals <- function(variable, deps = list(), priormarg = TRUE, indent = 0, verbose = TRUE, ...) { # deps necessary for ops functionality
	if (verbose) cat(rep("-", indent), "Checking", variable@name, "marginals", "...")
	varmar <- colnames(variable@data)[
		!grepl(paste("^", variable@name, "", sep=""), colnames(variable@data))&
		!colnames(variable@data) %in% c("Result", "Unit")
	]
	# all locs under observation/parameter index should be excluded
	varmar <- c(varmar, paste(variable@name, "Source", sep = "")) # Source is added 
	# by EvalOutput so it should always be in the initial list. 
	novarmar <- colnames(variable@data)[!colnames(variable@data) %in% varmar]
	if (priormarg & length(variable@marginal) > 0) {
		varmar <- unique(c(varmar, colnames(variable@output)[variable@marginal]))
		novarmar <- unique(c(novarmar, colnames(variable@output)[!colnames(variable@output) %in% varmar]))
	}
	if (length(deps) > 0) {
		for (i in deps) {
			if (class(i) == "ovariable") {
				varmar <- unique(c(varmar, colnames(i@output)[i@marginal]))
				novarmar <- unique(c(novarmar, colnames(i@output)[!i@marginal]))
			}
		}
	} else {
		for (i in as.character(variable@dependencies$Name)){
			if (class(get(i)) == "ovariable") {
				varmar <- unique(c(varmar, colnames(get(i)@output)[get(i)@marginal]))
				novarmar <- unique(c(novarmar, colnames(get(i)@output)[!get(i)@marginal]))
			}
		}
	}
	varmar <- varmar[!varmar %in% novarmar]
	variable@marginal <- colnames(variable@output) %in% c(varmar, "Iter")
	if (sum(variable@marginal) > 0) {
		if (verbose) cat("", paste(colnames(variable@output)[variable@marginal], collapse = ", "), "recognized as marginal(s).\n")
	} else {if (verbose) cat(" none recognized.\n")}
	#cat("done!\n")
	return(variable)
}