# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009
# Version 1.0
# Licence GPL v3

setGeneric("transition", function(x, transitionFunction, directions, ...) standardGeneric("transition"))

setMethod("transition", signature(x = "RasterLayer"), def = function(x, transitionFunction, directions, symm=TRUE, intervalBreaks=NULL)
		{
			if(class(transitionFunction)=="character") 
			{
				if(transitionFunction != "barriers" & transitionFunction != "areas")
				{
					stop("argument transitionFunction invalid")
				}
				if(transitionFunction=="barriers")
				{
					return(.barriers(x, directions, symm, intervalBreaks))
				}
				if(transitionFunction=="areas")
				{
					return(.areas(x, directions))
				}
			} else {
				return(.TfromR(x, transitionFunction, directions, symm))
			}
		}
)

.TfromR <- function(x, transitionFunction, directions, symm)
{
	tr <- new("TransitionLayer",
		nrows=as.integer(nrow(x)),
		ncols=as.integer(ncol(x)),
		extent=extent(x),
		crs=projection(x, asText=FALSE),
		transitionMatrix = Matrix(0,ncell(x),ncell(x)),
		transitionCells = 1:ncell(x))
	transitionMatr <- transitionMatrix(tr)
	Cells <- which(!is.na(getValues(x)))
	adj <- adjacent(x, cells=Cells, pairs=TRUE, target=Cells, directions=directions)
	if(symm){adj <- adj[adj[,1] < adj[,2],]}
	dataVals <- cbind(getValues(x)[adj[,1]],getValues(x)[adj[,2]])
	transition.values <- apply(dataVals,1,transitionFunction)
	if(!all(transition.values>=0)){warning("transition function gives negative values")}
	transitionMatr[adj] <- as.vector(transition.values)
	if(symm)
	{
		transitionMatr <- forceSymmetric(transitionMatr)
	}
	transitionMatrix(tr) <- transitionMatr
	matrixValues(tr) <- "conductance"
	return(tr)
}

.barriers <- function(x, directions, symm, intervalBreaks) {
	Xlayer <- new("TransitionLayer",
		nrows=as.integer(nrow(x)),
		ncols=as.integer(ncol(x)),
		extent=extent(x),
		crs=projection(x, asText=FALSE),
		transitionMatrix = Matrix(0,ncell(x),ncell(x)),
		transitionCells = 1:ncell(x))
	matrixValues(Xlayer) <- "resistance"
	Xstack <- as(Xlayer, "TransitionStack") * 0
	#Xstack@transition <- vector(list,...)
	
	if(x@data@isfactor) {

		vals <- unlist(x@data@attributes[[1]])
		n <- length(vals)
		
		if(symm)
		{
			maxn <- (n^2 - n)/2
			for(i in 1:maxn)
			{
				j <- .matrIndex(i,n)
				XlayerNew <- Xlayer
				cells1 <- which(getValues(x) == vals[j[1]])
				cells2 <- which(getValues(x) == vals[j[2]])
				adj1 <- adjacent(x, cells=cells1, pairs=TRUE, target=cells2, directions=directions)
				adj2 <- adjacent(x, cells=cells2, pairs=TRUE, target=cells1, directions=directions)
				adj <- rbind(adj1,adj2)
				XlayerNew[adj] <- 1
				Xstack <- stack(Xstack, XlayerNew)
			}
		} else {
			maxn <- (n^2 - n)/2
			for(i in 1:maxn)
			{
				j <- .matrIndex(i,n)
				XlayerNew1 <- Xlayer
				XlayerNew2 <- Xlayer
				cells1 <- which(getValues(x) == vals[j[1]])
				cells2 <- which(getValues(x) == vals[j[2]])
				adj1 <- adjacent(x, cells=cells1, pairs=TRUE, target=cells2, directions=directions)
				adj2 <- adjacent(x, cells=cells2, pairs=TRUE, target=cells1, directions=directions)
				XlayerNew1[adj1] <- 1
				XlayerNew2[adj2] <- 1				
				Xstack <- stack(Xstack, XlayerNew1, XlayerNew2)
			}
		}
	
	} else {
	
		Xmin <- transition(x, min, directions)
		Xmax <- transition(x, max, directions)
		index1 <- adjacent(x, cells=1:ncell(x), pairs=TRUE, target=1:ncell(x), directions=directions)
		XminVals <- Xmin[index1]
		XmaxVals <- Xmax[index1]

		if(symm == TRUE)
		{
			for(i in 1:length(intervalBreaks))
			{
				index2 <- index1[XminVals < intervalBreaks[i] & XmaxVals > intervalBreaks[i],]
				XlayerNew <- Xlayer
				XlayerNew[index2] <- 1
				Xstack <- stack(Xstack,XlayerNew)
			}
		}
		if(symm=="up" | symm=="down"){stop("not implemented yet")}
	}
	
	Xstack <- Xstack[[2:nlayers(Xstack)]]	
	return(Xstack)
}


.areas <- function(x, directions) {

	Xlayer <- new("TransitionLayer",
		nrows=as.integer(nrow(x)),
		ncols=as.integer(ncol(x)),
		extent=extent(x),
		crs=projection(x, asText=FALSE),
		transitionMatrix = Matrix(0,ncell(x),ncell(x)),
		transitionCells = 1:ncell(x))
	matrixValues(Xlayer) <- "resistance"
	Xstack <- as(Xlayer, "TransitionStack") * 0
	#Xstack@transition <- vector(list,...)
	
	if(x@data@isfactor) {

		vals <- unlist(x@data@attributes[[1]])
		n <- length(vals)
		
			for(i in 1:n)
			{
				transitionFunction <- function(v) {return(sum(v == i) / 2)}
				XlayerNew <- .TfromR(x, transitionFunction, directions, symm=TRUE)
				Xstack <- stack(Xstack,XlayerNew)
			}
			
	} else {
		warning("not yet implemented for raster with non-factor variables. Contact author.")
	}
	Xstack <- Xstack[[2:nlayers(Xstack)]]	
	return(Xstack)
}



setMethod("transition", signature(x = "RasterBrick"), def = function(x, transitionFunction="mahal", directions)
		{
			if(transitionFunction != "mahal")
			{
				stop("only Mahalanobis distance method implemented for RasterBrick")
			}
			xy <- cbind(1:ncell(x),getValues(x))
			xy <- na.omit(xy)
			dataCells <- xy[,1]
			adj <- adjacent(x, cells=dataCells, pairs=TRUE, target=dataCells, directions=directions)
			x.minus.y <- xy[adj[,1],-1]-xy[adj[,2],-1]
			cov.inv <- solve(cov(xy[,-1]))
			mahaldistance <- apply(x.minus.y,1,function(x){sqrt((x%*%cov.inv)%*%x)})
			mahaldistance <- mean(mahaldistance)/(mahaldistance+mean(mahaldistance))
			transitiondsC <- new("dsCMatrix", 
					p = as.integer(rep(0,ncell(x)+1)),
					Dim = as.integer(c(ncell(x),ncell(x))),
					Dimnames = list(as.character(1:ncell(x)),as.character(1:ncell(x)))
			)
			transitiondsC[adj] <- mahaldistance
			tr <- new("TransitionLayer",nrows=as.integer(nrow(x)),ncols=
			as.integer(ncol(x)),extent = extent(x),
			crs=projection(x, asText=FALSE), matrixValues="conductance", transitionMatrix = transitiondsC)
			return(tr)
		}
)