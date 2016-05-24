# Author: Jacob van Etten jacobvanetten@yahoo.com
# Bioversity International

setGeneric("transitionMatrix", function(x, inflate) standardGeneric("transitionMatrix"))

setMethod ("transitionMatrix", signature(x = "TransitionLayer", inflate="missing"),
	function(x)
	{
		transitionMatrix(x=x, inflate=TRUE)
	}
)

setMethod ("transitionMatrix", signature(x = "TransitionLayer", inflate="logical"),
	function(x, inflate)
	{
		.tr(x, inflate, ncell(x))
	}
)

setMethod ("transitionMatrix", signature(x = "TransitionData", inflate="missing"),
	function(x)
	{
		.tr(x=x, inflate=FALSE, nc=0)
	}
)

.tr <- function(x, inflate, nc)
{
	if(inflate & length(transitionCells(x)) != ncell(x))
	{
		tr <- Matrix(0, nc,nc)
		cells <- transitionCells(x)
		tr[cells,cells] <- x@transitionMatrix
	}
	if(!inflate | length(transitionCells(x)) == nc)
	{
		tr <- x@transitionMatrix
	}
	return(tr)
}


setGeneric("transitionMatrix<-", function(x, value) standardGeneric("transitionMatrix<-"))

setReplaceMethod ("transitionMatrix", signature(x = "TransitionLayer", value = "sparseMatrix"),
	function(x, value)
	{
		if(dim(value)[1] != dim(value)[2]){stop("sparse matrix has to be square")}
		if(dim(value)[1] != ncell(x)[2]){stop("sparse matrix has to have ncell(x) rows and columns")}
		x@transitionMatrix <- value
		x@transitionCells <- 1:ncell(x)
		return(x)
	}
)

setGeneric("transitionCells", function(x) standardGeneric("transitionCells"))

setMethod ("transitionCells", signature(x = "TransitionLayer"),
	function(x)
	{
		return(x@transitionCells)
	}
)

setMethod ("transitionCells", signature(x = "TransitionData"),
	function(x)
	{
		return(x@transitionCells)
	}
)

setGeneric("matrixValues", function(x) standardGeneric("matrixValues"))

setMethod ("matrixValues", signature(x = "TransitionLayer"),
	function(x){x@matrixValues}
)

setMethod ("matrixValues", signature(x = "TransitionStack"),
	function(x){stop("not implemented yet")}
)

setGeneric("matrixValues<-", function(x, value) standardGeneric("matrixValues<-"))

setReplaceMethod ("matrixValues", signature(x = "TransitionLayer", value = "character"),
	function(x, value){
		if (value == "resistance" | value == "conductance") 
		{
			x@matrixValues <- value
			return(x)
		}
		else {stop("matrixValues can only be set to resistance or conductance")}
	}
)

setMethod("[", signature(x = "TransitionLayer", i="numeric", j="numeric", drop="missing"), function(x,i,j)
	{
		tm <- transitionMatrix(x)
		tm <- tm[i,j]
		return(tm)
	}
)

setMethod("[", signature(x = "TransitionLayer", i="matrix", j="missing", drop="missing"), function(x,i)
	{
		tm <- transitionMatrix(x)
		tm <- tm[i]
		return(tm)
	}
)

setMethod("[<-", signature(x = "TransitionLayer", i="matrix", j="missing", value="ANY"),
		function(x, i, value){
			tm <- transitionMatrix(x)
			tm[i] <- value
			x@transitionMatrix <- tm
			return(x)
		}
)

setMethod("[<-", signature(x = "TransitionLayer", i="numeric", j="numeric", value="ANY"),
		function(x, i, j, value)
		{
			tm <- transitionMatrix(x)
			tm[i,j] <- value
			transitionMatrix(x) <- tm
			return(x)
		}
)

setGeneric("transitionMatrix<-", function(x, value) standardGeneric("transitionMatrix<-"))

setReplaceMethod ("transitionMatrix", signature(x = "TransitionLayer", value = "sparseMatrix"),
	function(x, value){
		if(dim(value)[1] != dim(value)[2]){stop("sparse matrix has to be square")}
		if(dim(value)[1] == ncell(x)){x@transitionMatrix <- value}
		else
		{
			if(dim(value)[1] == length(transitionCells(x)))
			{
				trC <- transitionCells(x)
				tr <- Matrix(0,ncell(x),ncell(x))
				tr[trC,trC] <- value
				x@transitionMatrix <- tr
			}
			else{stop("value is of wrong dimensions; either ncell(transition) or length(transitionCells(transition))")}
		}
		return(x)
	}
)

setMethod('nlayers', signature(x='TransitionStack'), 
	function(x)
	{
		return(length(x@transition)) 
    }
)


setMethod("[[", signature(x = "TransitionStack", i="numeric", j="missing"), function(x,i)
	{
		if (!(all(i %in% 1:nlayers(x)))){stop("indices should correspond to layers")}
		else
		{
			if(length(i)==1)
			{
				result <- new("TransitionLayer", nrows=as.integer(nrow(x)),ncols = as.integer(ncol(x)), extent = extent(c(xmin(x), xmax(x),
				ymin(x), ymax(x))), crs=projection(x, asText=FALSE))
				result@transitionMatrix <- x@transition[[i]]@transitionMatrix
				result@transitionCells <- x@transition[[i]]@transitionCells
				result@matrixValues <- x@transition[[i]]@matrixValues
			}			
			if(length(i)>1)
			{
				result <- x
				result@transition <- x@transition[i]
			}
		}
		return(result)
	}
)

setMethod("[[<-", signature(x = "TransitionStack", i="numeric", j="missing", value="TransitionData"), function(x,i, value)
	{
		x@transition[[i]] <- value
		return(x)
	}
)

setGeneric("transitionData", function(x) standardGeneric("transitionData"))

setMethod ("transitionData", signature(x = "TransitionLayer"),
	function(x){
		as(x, "TransitionData")
	}
)

setMethod ("transitionData", signature(x = "TransitionStack"),
	function(x){
		x@transition
	}
)