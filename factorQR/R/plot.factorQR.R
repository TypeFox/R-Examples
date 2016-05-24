plot.factorQR <- function(x, ...){
	dots = list(...)

	if(is.null(dots$whichPlot)) dots$whichPlot <- "LambdaQ"
	if(is.null(dots$plotBurn)) dots$plotBurn <- 0

		param <- x$param
		n.param <- length(param[1,])
		nFact <- x$nFact
		nSamp <- x$nSamp
		burn <- x$burn
		thin <- x$thin
		nFactorX <- x$nFactorX
		betLen <- x$betLen
		nStored <- dim(param)[1]
		whichPlot <- dots$whichPlot
		nFactInt <- x$nFactInt

	if(is.numeric(whichPlot)){
		sampledVals <- param[, whichPlot]
		nToPlot <- length(whichPlot)
		whichNames <- whichPlot
		}
	else if (whichPlot == "LambdaQ"){
		sampledVals <- param[,1:nFactInt]
		nToPlot <- nFactInt
		whichNames <- 1:nFactInt
		}
	else if (whichPlot == "Lambda"){
		sampledVals <- param[, 1:(nFactInt + nFactorX)]
		nToPlot <- nFact + nFactorX
		whichNames <- 1:(nFactInt + nFactorX)}
	else if (whichPlot == "beta"){
		sampledVals <- param[, (nFactInt + nFactorX + 1) : (nFactInt + nFactorX + betLen)]
		nToPlot <- length((nFactInt + nFactorX + 1) : (nFactInt + nFactorX + betLen))
		whichNames <- (nFactInt + nFactorX + 1) : (nFactInt + nFactorX + betLen)}
	else if (whichPlot == "variances"){
		sampledVals <- param[, -((nFactInt + nFactorX + 1) : (nFactInt + nFactorX + betLen))]
		nToPlot <- n.param - (nFactInt + nFactorX + betLen)
		whichNames <- (nFactInt + nFactorX + betLen+1):n.param}
	else {
	stop("`whichPlot' must be NULL, `LambdaQ', `Lambda', `beta', `variances' or a numeric vector of the columns of param to plot.")}
	
		
		sampledVals <- matrix(sampledVals, nc=1)
		comp <- rep(colnames(param)[whichNames], each=nStored)
		iteration <- rep(seq((burn+1), nSamp, length=nStored), nToPlot)
		ptData = data.frame(sampledVals, comp, iteration)
		pt <- xyplot(sampledVals ~ iteration|comp, data=ptData, xlab= if(thin>0) "Thinned Iteration" else "Iteration", ylab="Sampled Values", strip=FALSE, strip.left=TRUE, layout = c(1,nToPlot), type=c("g","l"), subset=iteration >= dots$plotBurn, scales=list(y=list(relation="free")))
		return(pt)
		
	}
