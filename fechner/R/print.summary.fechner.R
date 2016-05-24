#####################################################################################
##print method for objects of class "summary.fechner"; NOTE: this function is used###
##for printing summary information about objects of class "fechner" computed by######
##function summary.fechner; it is called through using the generic function summary##
#####################################################################################
print.summary.fechner <-
function(x, ...){
	# x: an object of class "summary.fechner", obtained from a call to the internal function summary.fechner
	#    through using the generic function summary
	# ...: further arguments to be passed to or from other methods; they are ignored in this function

	cat("\n")
	cat("number of stimuli pairs used for comparison:", length(x$pairs.used.for.comparison[, 1]), "\n")
	cat("\n")
	cat("summary of corresponding S-index values:\n")
	print(summary(x$pairs.used.for.comparison[, 2]))
	cat("\n")
	cat("summary of corresponding Fechnerian distance G values:\n")
	print(summary(x$pairs.used.for.comparison[, 3]))
	cat("\n")
	cat("Pearson correlation:", x$Pearson.correlation, "\n")
	cat("\n")
	cat("C-index:", x$C.index, "\n")
	cat("\n")
	cat("comparison level:", x$comparison.level, "\n")
	cat("\n")
	invisible(x)
}
