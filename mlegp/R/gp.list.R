`gp.list` <-
function(..., param.names = NULL, UD = NULL, gp.names = NULL) {
		
	if (is.list(...)) {
		x <- as.list(...)
	}
	else {
		x <- list(...)
		if (length(x) == 1) {
			stop("Cannot create list for 1 gp")
		}
	}
	numDims = lapply(x, gpDim)
	if (!all(unlist(lapply(numDims, "==", numDims[[1]])))) {
		stop("Different number of dimensions in at least one gp")
	}
	numObs = lapply(x, gpObs)	
	if (!all(unlist(lapply(numObs, "==", numObs[[1]])))) {
		stop("Different number of observations in at least one gp")
	}
	
	numGPs = length(x)

	if (is.null(param.names)) {
		param.names = x[[1]]$params
	}

	if (is.null(gp.names)) {
		gp.names = paste("gp #", 1:numGPs)
	}

	if (length(gp.names) != numGPs) {
		stop("length of gp.names must match number of GPs")
	}

	x$names = gp.names
	x$params = param.names
	x$numGPs = numGPs
	x$numDim = numDims[[1]]
	x$numObs = numObs[[1]]
	x$UD = UD

	attr(x, "class") <- "gp.list"
        x
}

