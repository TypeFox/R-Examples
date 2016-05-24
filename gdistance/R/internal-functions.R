# Author: Jacob van Etten, jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009
# Version 1.0
# Licence GPL v3

.coordsToMatrix <- function(Coords)
{
	if(class(Coords) == "numeric")
	{
		if(length(Coords) == 2) {Coords <- t(as.matrix(Coords))} 
		else{stop("coordinates given as a vector, but the vector does not have a length of two")}
	}
	
	if(class(Coords) == "matrix")
	{
		if(!(ncol(Coords) == 2)){stop("coordinates given as a matrix, but the matrix does not have two columns")}
	}	

	if(inherits(Coords, "SpatialPoints"))  
	{
		Coords <- coordinates(Coords)
	}
	return(Coords)
}

.connected.components <- function(x)
{
	adj.graph <- graph.adjacency(transitionMatrix(x))
	clustermembership <- cbind(1:ncell(x),as.integer(clusters(adj.graph)$membership)+1)
	return(clustermembership)
}

.current <- function(L, Lr, A, n, indexFrom, indexTo) 
{
	C <- 1e-300 * n #This should avoid too big floating points as "Voltage differences"
	e <- matrix(0, ncol=1, nrow=n)
	e[indexFrom,] <- C
 	e[indexTo,] <- -C
	x <- solve(Lr,e)
	x <- as.vector(x)
	Lplusallrows <- c(x,x[length(x)]) / C
	V1 <- A * Lplusallrows
	V2 <- t(t(A) * Lplusallrows)
	V <- abs(V1 - V2)
	Current <- colSums(V * -L)/2 #I = V * Conductance
	Current[indexFrom] <- 1
	Current[indexTo] <- 1
	return(Current)
}

.currentR <- function(L, Lr, A, n, indexFrom, indexTo)
{
	lf <- length(indexFrom)
	lt <- length(indexTo)
	C <- 1e-300 * n
	Cf <- C / lf #This should avoid too big floating points as "Voltage differences"
	Ct <- C / lt
	e <- matrix(0, ncol=1, nrow=n)
	e[indexFrom,] <- Cf
 	e[indexTo,] <- -Ct
	x <- solve(Lr,e)
	x <- as.vector(x)
	Lplusallrows <- c(x,x[length(x)]) / C
	V1 <- A * Lplusallrows
	V2 <- t(t(A) * Lplusallrows)
	V <- abs(V1 - V2)
	Current <- colSums(V * -L)/2 #I = V * Conductance
	Current[indexFrom] <- 1
	Current[indexTo] <- 1
	return(Current)
}

.potential <- function(L, Lr, A, n, indexFrom, indexTo) 
{
	C <- 1e-300 * n #This should avoid too big floating points as "Voltage differences"
	e <- matrix(0, ncol=1, nrow=n)
	e[indexFrom,] <- C
 	e[indexTo,] <- -C
	x <- solve(Lr,e)
	x <- as.vector(x)
	Lplusallrows <- c(x,x[length(x)]) / C
	V1 <- A * Lplusallrows
	V2 <- t(t(A) * Lplusallrows)
	V <- abs(V1 - V2)
	return(V)
}

.currentM <- function(L, Lr, A, n, indexFrom, indexTo, index) 
{
	C <- 1e-300 * n #This should avoid too big floating points as "Voltage differences"
	e <- matrix(0, ncol=1, nrow=n)
	e[indexFrom,] <- C
 	e[indexTo,] <- -C
	x <- solve(Lr,e)
	x <- as.vector(x)
	Lplusallrows <- c(x,x[length(x)]) / C
	V1 <- A * Lplusallrows
	V2 <- t(t(A) * Lplusallrows)
	V <- abs(V1 - V2)
	Current <- V[index] * -L[index] #I = V * Conductance
	return(Current)
}

.Laplacian <- function(x) 
{
	Laplacian <- Diagonal(x = colSums(transitionMatrix(x, inflate=FALSE))) - transitionMatrix(x, inflate=FALSE)
	Laplacian <- as(Laplacian, "symmetricMatrix")
	return(Laplacian)
}

.transitionSolidify <- function(x)
{
	selection <- which(rowMeans(transitionMatrix(x,inflate=FALSE))>1e-300)
	x@transitionCells <- x@transitionCells[selection]
	x@transitionMatrix <- transitionMatrix(x,inflate=FALSE)[selection,selection]
	return(x)
}

#determine place in dist vector given place in dist matrix -- from gdistanalyst
.distIndex <- function(i,j,n){n*(j-1) - j*(j-1)/2 + i-j}

#determine place in dist matrix given place in dist vector -- from gdistanalyst -- should be possible speed up!
.matrIndex <- function(i,n){
	cc <- cumsum(seq((n-1),1))
	out <- matrix(nrow=length(i),ncol=2)
	for(index in 1:length(i))
	{
		out[index,2] <- min(which((cc-i[index])>=0))
		out[index,1] <- -c(0,cc)[out[index,2]] + i[index] + out[index,2]
	}
	return(out)
}


