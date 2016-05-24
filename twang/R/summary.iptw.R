summary.iptw <- function (object, ...){
	nFits <- object$nFits
	summaryList <- vector(mode = "list", length = nFits)
	for (i in 1:nFits) {
		summaryList[[i]] <- summary(object$psList[[i]], ...)
	}

	retObj <- list(summaryList = summaryList, nFit = object$nFit, uniqueTimes = object$uniqueTimes)

	class(retObj) <- "summary.iptw"
	return(retObj)
}

