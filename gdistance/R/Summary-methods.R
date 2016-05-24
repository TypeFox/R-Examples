# Author: Jacob van Etten
# Date: April 2011
# Version 1.1.2
# Licence GPL v3

setMethod("Summary", signature(x='TransitionStack'),
	function(x, ..., na.rm=FALSE){
		objectList <- list(x, ...)
		if(!is.null(objectList))
		{
			if(length(objectList)>1) warning("operations with more than two Transition* objects not implemented; use stack()")
		}
		call <- sys.call()
        fun <- as.character(call[[1L]])
		result <- x[[1]] * 0
		trResult <- .rowWiseFun(x,result,fun)
		transitionMatrix(result) <- trResult
		return(result)
	}
)

.rowWiseFun <- function(x, result, fun)
{
		trResult <- transitionMatrix(result)
		for(i in 1: nrow(trResult))
		{
			rows <- x[[1]]@transitionMatrix[i,]
			for(j in 2:nlayers(x)) {rows <- rbind(rows, x[[j]]@transitionMatrix[i,])}
			trResult[i,] <- apply(rows, 2, fun)
		}
		return(trResult)
}

setMethod("Summary", signature(x = "TransitionLayer"),
		function(x, ..., na.rm)
		{
			objectList <- list(x, ...)
			if(!is.null(objectList))
			{
				if(length(objectList)>1)
				{
					stop("operations with more than two Transition* objects not implemented; use stack()")
				}
				else{result <- callGeneric(x@transitionMatrix, objectList[[1]]@transitionMatrix, na.rm=na.rm)}
			}
			else{result <- callGeneric(x@transitionMatrix, na.rm=na.rm)}

			return(result)
		}
)

setMethod("sum", signature(x='TransitionStack'),
	function(x, ..., na.rm=FALSE){
		objectList <- list(x, ...)
		if(!is.null(objectList))
		{
			if(length(objectList)>1) warning("operations with more than two Transition* objects not implemented; use stack()")
		}
		result <- x[[1]] * 0
		trResult <- .MatrixSum(x,result) 
		transitionMatrix(result) <- trResult
		return(result)
	}
)

.rowWiseSum <- function(x, result)
{
		trResult <- transitionMatrix(result)
		for(i in 1: nrow(trResult))
		{
			rows <- x[[1]]@transitionMatrix[i,]
			for(j in 2:nlayers(x)) {rows <- rbind(rows, x[[j]]@transitionMatrix[i,])}
			trResult[i,] <- colSums(rows)
		}
		return(trResult)
}

.MatrixSum <- function(x, result)
{
		trResult <- transitionMatrix(result)
		for(i in 1: nlayers(x))
		{
			m <- x[[i]]@transitionMatrix
			trResult <- trResult + m
		}
		return(trResult)
}

setMethod("mean", signature(x='TransitionStack'),
	function(x, ..., na.rm=FALSE){
		objectList <- list(x, ...)
		if(!is.null(objectList))
		{
			if(length(objectList)>1) warning("operations with more than two Transition* objects not implemented; use stack()")
		}
		result <- x[[1]] * 0
		trResult <- .rowWiseMean(x,result)
		transitionMatrix(result) <- trResult
		return(result)
	}
)

.rowWiseMean <- function(x, result)
{
		trResult <- transitionMatrix(result)
		for(i in 1: nrow(trResult))
		{
			rows <- x[[1]]@transitionMatrix[i,]
			for(j in 2:nlayers(x)) {rows <- rbind(rows, x[[j]]@transitionMatrix[i,])}
			trResult[i,] <- colMeans(rows)
		}
		return(trResult)
}