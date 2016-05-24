ida <- function(y, yp, times, package, resfunc, jacfunc = NULL, rootfunc = NULL, fndata = NULL, numroots = 0, rtol = 1e-6, 
				atol = 1e-6, maxnumsteps = 500, maxstep = 0, verbose = FALSE, lasttime = FALSE) {
	
	if (!is.numeric(y)) stop("Error: 'y' must be numeric")
	if (!is.numeric(yp)) stop("Error: 'yp' must be numeric")
	if (!is.numeric(times)) stop("Error: 'times' must be numeric")
	if (!is.character(package)) stop("Error: 'package' must be a character vector")
	if (!is.character(resfunc)) stop("Error: 'resfunc' must be a character vector")
	if (!is.null(jacfunc) && !is.character(jacfunc)) stop("Error: 'jacfunc' must be a character vector")
	if (!is.null(rootfunc) && !is.character(rootfunc)) stop("Error: 'rootfunc' must be a character vector")
	if (!is.null(rootfunc) && numroots <= 0) stop("Error: numroots must be greater than 0");
	if (!is.numeric(numroots)) stop("Error: 'numroots' must be numeric")
	if (!is.numeric(rtol)) stop("Error: 'rtol' must be numeric")
	if (!is.numeric(atol)) stop("Error: 'atol' must be numeric")
	if (!is.numeric(maxnumsteps)) stop("Error: 'maxsteps' must be numeric")
	if (!is.numeric(maxstep)) stop("Error: 'maxsteps' must be numeric")
	if (!is.null(fndata) && !is.numeric(fndata)) stop("Error: Data arguments must be numeric");
	
	refunc = getNativeSymbolInfo(resfunc,PACKAGE=package)$address
	jfunc = NULL
	if (!is.null(jacfunc)) jfunc = getNativeSymbolInfo(jacfunc, PACKAGE=package)$address
	rofunc = NULL
	if (!is.null(rootfunc)) rofunc = getNativeSymbolInfo(rootfunc, PACKAGE=package)$address
	solutions = c()
	solutions = .Call("ida",PACKAGE="Rsundials", as.double(y), as.double(yp), as.double(times), refunc, jfunc, rofunc, as.double(fndata),
					as.integer(numroots), as.double(rtol), as.double(atol), as.integer(maxnumsteps), as.integer(maxstep), as.integer(verbose), as.integer(lasttime))
	
	if (lasttime == TRUE) rows = 1
	else rows = length(times)
	
	solutions <- matrix(solutions, rows)
	col <- c()
	for (i in 1:length(y)) col[i] <- paste("y",i,sep="")
	if (lasttime == TRUE) dimnames(solutions) <- list(times[length(times)],col)
	else dimnames(solutions) <- list(times,col)
	
	solutions
}
