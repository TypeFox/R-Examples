# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute, Bioversity International
# Date :  January 2009, normalization added May 2012
# Version beta
# Licence GPL v3

#TODO check if coordinate systems are equal.
#TODO check if bounding box of coordinates falls inside bb of transition
#TODO coordinates in same cell: distance = 0

setGeneric("commuteDistance", function(x, coords) standardGeneric("commuteDistance"))

setMethod("commuteDistance", signature(x = "TransitionLayer", coords = "Coords"), def = function(x, coords) 
	{
		return(.rD(x, coords))
	}
)

setMethod("commuteDistance", signature(x = "TransitionLayer", coords = "Coords"), def = function(x, coords) 
{
  return(.rD(x, coords))
}
)

.rD <- function(x, coords){
		if(class(transitionMatrix(x)) != "dsCMatrix"){stop("symmetric transition matrix required (dsCMatrix) in TransitionLayer object x")}
		coords <- .coordsToMatrix(coords)

		rd <- matrix(NA,nrow=length(coords[,1]),ncol=length(coords[,1]))
		rownames(rd) <- rownames(coords)
		colnames(rd) <- rownames(coords)
		allFromCells <- cellFromXY(x, coords)
		
		if(!all(!is.na(allFromCells))){
			warning("some coordinates not found and omitted")
			allFromCells <- allFromCells[!is.na(allFromCells)]
		}

		x <- .transitionSolidify(x)		
		fromCells <- allFromCells[allFromCells %in% transitionCells(x)]
		if (length(fromCells) < length(allFromCells)) 
		{
			warning(length(fromCells)," out of ",length(allFromCells)," locations were found in the fully connected transition matrix. NAs introduced.")
		}
		else{}
		fromCells <- unique(allFromCells)

  
	Lr <- .Laplacian(x)
	n <- max(Lr@Dim)
	Lr <- Lr[-n,-n]
	C <- 1e-300 * n #This should avoid too big floating points as "Voltage differences", but give a number that can still be divided by n
	Lplus <- matrix(ncol=length(fromCells),nrow=length(fromCells))
	index <- match(fromCells,transitionCells(x))
	#Lr <- Cholesky(Lr)
	for (i in 1:length(fromCells))
	{
		ei <- matrix(-C/n, ncol=1, nrow=n-1)
		ei[index[i],] <- C-(C/n)
		xi <- solve(Lr,ei)
		xi <- as.vector(xi)
		Lplusallrows <- c(xi-sum(xi/n),(sum(xi)/n)) 
		Lplus[,i] <- as.vector(Lplusallrows)[index]
	}
	Lplus <- Lplus / C
	rdSS <- (-2*Lplus + matrix(diag(Lplus),nrow=length(fromCells),ncol=length(fromCells)) 
			+ t(matrix(diag(Lplus),nrow=length(fromCells),ncol=length(fromCells)))) 

    Volume <- sum(transitionMatrix(x))
    rdSS <- rdSS * Volume

    index1 <- which(allFromCells %in% fromCells)
		index2 <- match(allFromCells[allFromCells %in% fromCells],fromCells)
		rd[index1,index1] <- rdSS[index2,index2]
		rd <- as.dist(rd)
		attr(rd, "method") <- "commute"
		return(rd)
}