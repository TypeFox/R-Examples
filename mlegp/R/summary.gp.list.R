`summary.gp.list` <-
function(object, nums = NULL, ...) {

	if (is.null(nums)) {
		cat ("\nnum GPs: ")
		cat (object$numGPs)
		cat ("\nTotal observations (per GP): ")
		cat (object$numObs)
		cat ("\nDimensions: ")
		cat (object$numDim)
		cat ("\n\n")
		return 
	}

	for (i in nums) {
		cat ("===== ")	
		cat (object$names[i])
		cat (" =====")
		cat (" \n")
		summary(object[[i]])
		cat("\n")
	}
}

