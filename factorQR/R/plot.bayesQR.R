plot.bayesQR <- function(x, ...){
	dots = list(...)

	if(is.null(dots$whichPlot)) dots$whichPlot <- "beta"
	if(is.null(dots$plotBurn)) dots$plotBurn <- 0

		param <- x$param
		n.param <- length(param[1,])
		nSamp <- x$nSamp
		burn <- x$burn
		thin <- x$thin
		betLen <- x$betLen
		nStored <- dim(param)[1]
		whichPlot <- dots$whichPlot

	if(is.numeric(whichPlot)){
		sampledVals <- param[, whichPlot]
		nToPlot <- length(whichPlot)
		whichNames <- whichPlot
		}
	else if (whichPlot == "tau"){
		sampledVals <- param[,n.param]
		nToPlot <- 1
		whichNames <- n.param
		}
	else if (whichPlot == "beta"){
		sampledVals <- param[, -n.param]
		nToPlot <- n.param - 1
		whichNames <- 1:(n.param - 1)}
	else {
	stop("`whichPlot' must be NULL, `beta', `tau' or a numeric vector of the columns of `param' to plot.")}
	
		
		sampledVals <- matrix(sampledVals, nc=1)
		comp <- rep(colnames(param)[whichNames], each=nStored)
		iteration <- rep(seq((burn+1), nSamp, length=nStored), nToPlot)
		ptData = data.frame(sampledVals, comp, iteration)
		pt <- xyplot(sampledVals ~ iteration|comp, data=ptData, xlab= if(thin>0) "Thinned Iteration" else "Iteration", ylab="Sampled Values", strip=FALSE, strip.left=TRUE, layout = c(1,nToPlot), type=c("g","l"), subset=iteration >= dots$plotBurn, scales=list(y=list(relation="free")))
		return(pt)
		
	}
