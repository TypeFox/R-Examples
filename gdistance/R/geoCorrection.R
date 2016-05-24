# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009
# Version beta
# Licence GPL v3


setGeneric("geoCorrection", function(x, type, ...) standardGeneric("geoCorrection"))

setMethod("geoCorrection", signature(x = "TransitionLayer", type="missing"), def = function(x, multpl=FALSE, scl=FALSE)
	{
		return(geoCorrection(x, type="c", multpl, scl))
	}
)


setMethod("geoCorrection", signature(x = "TransitionLayer", type="character"), def = function(x, type, multpl=FALSE, scl=FALSE)
	{
		scaleValue <- 1
		if(scl)
		{
			midpoint <- c(mean(c(xmin(x),xmax(x))),mean(c(ymin(x),ymax(x))))
			scaleValue <- pointDistance(midpoint,midpoint+c(xres(x),0), longlat=isLonLat(x)) 
		}
		if(isLonLat(x))
		{
			if (type != "c" & type != "r"){stop("type can only be c or r")}
			#if (type == "r" & matrixValues(x) != "conductance"){stop("matrix of Transition object must have conductance values")}
			adj <- adjacencyFromTransition(x)
			correction <- cbind(xyFromCell(x,adj[,1]),xyFromCell(x,adj[,2]))
			if(matrixValues(x) == "conductance") {correctionValues <- 1/(pointDistance(correction[,1:2],correction[,3:4],longlat=TRUE)/scaleValue)}
			if(matrixValues(x) == "resistance") {correctionValues <- pointDistance(correction[,1:2],correction[,3:4],longlat=TRUE)/scaleValue}
			if (type=="r")
			{
				rows <- rowFromCell(x,adj[,1]) != rowFromCell(x,adj[,2])
				if(matrixValues(x) == "conductance") {corrFactor <- cos((pi/180) * rowMeans(cbind(correction[rows,2],correction[rows,4])))} #low near the poles
				if(matrixValues(x) == "resistance") {corrFactor <- 1 / (cos((pi/180) * rowMeans(cbind(correction[rows,2],correction[rows,4]))))} #high near the poles
				correctionValues[rows] <- correctionValues[rows] * corrFactor #makes conductance lower in N-S direction towards the poles
			}
		} else {
			adj <- adjacencyFromTransition(x)
			correction <- cbind(xyFromCell(x,adj[,1]),xyFromCell(x,adj[,2]))
			if(matrixValues(x) == "conductance") {correctionValues <- 1/(pointDistance(correction[,1:2],correction[,3:4],longlat=FALSE)/scaleValue)}
			if(matrixValues(x) == "resistance") {correctionValues <- pointDistance(correction[,1:2],correction[,3:4],longlat=FALSE)/scaleValue}
		}
		i <- as.integer(adj[,1] - 1)
		j <- as.integer(adj[,2] - 1)
		xv <- as.vector(correctionValues) #check for Inf values!
		dims <- ncell(x)
		correctionMatrix <- new("dgTMatrix", i = i, j = j, x = xv, Dim = as.integer(c(dims,dims)))
		correctionMatrix <- (as(correctionMatrix,"sparseMatrix"))
		if(class(transitionMatrix(x)) == "dsCMatrix"){correctionMatrix <- forceSymmetric(correctionMatrix)} #isSymmetric?
		if(!multpl) 
		{
			transitionCorrected <- correctionMatrix * transitionMatrix(x)
			transitionMatrix(x) <- transitionCorrected
			return(x)
		}	
		if(multpl)
		{
			transitionMatrix(x) <- correctionMatrix
			return(x)
		}
	}
)