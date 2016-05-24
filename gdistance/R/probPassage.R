# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009, code added January 2010
# Version 1.0
# Licence GPL v3

#check if Transition and RasterLayers coincide, etc.

setGeneric("passage", function(x, origin, goal, theta, ...) standardGeneric("passage"))

setMethod("passage", signature(x = "TransitionLayer", origin = "Coords", goal = "Coords", theta="missing"), def = function(x, origin, goal, totalNet="net", output="RasterLayer")
	{
		.checkInputsPassage(totalNet, output)
		
		origin <- .coordsToMatrix(origin)
		goal <- .coordsToMatrix(goal)
		
		if(totalNet=="net" & output=="RasterLayer")
		{
			x <- .transitionSolidify(x)
			tc <- transitionCells(x)
			cellnri <- cellFromXY(x, origin)
			cellnrj <- cellFromXY(x, goal)
			ci <- match(cellnri,tc)
			cj <- match(cellnrj,tc)
			result <- .flowMap(x, ci, cj, tc)
		}
		else{stop("no method available -- try a low value of theta instead")}
		return(result)
	}
)

setMethod("passage", signature(x = "TransitionLayer", origin = "RasterLayer", goal = "RasterLayer", theta="missing"), def = function(x, origin, goal, totalNet="net", output="RasterLayer")
	{
		.checkInputsPassage(totalNet, output)
		
		if(totalNet=="net" & output=="RasterLayer")
		{
			x <- .transitionSolidify(x)
			tc <- transitionCells(x)
			ci <- which(getValues(origin))
			cj <- which(getValues(goal))
			result <- .flowMap(x, ci, cj, tc)
		}
		else{stop("no method available -- try a low value of theta instead")}
		return(result)
	}
)

.checkInputsPassage <- function(totalNet, output)
{
	if(!(totalNet %in% c("total","net"))){stop("totalNet should be either total or net")}
	if(!(output %in% c("RasterLayer","TransitionLayer"))){stop("output should be either RasterLayer or TransitionLayer")}
} 

.flowMap <- function(x, indexOrigin, indexGoal, tc)
{
	L <- .Laplacian(x)
	Lr <- L[-dim(L)[1],-dim(L)[1]]
	A <- as(L,"lMatrix")
	A <- as(A,"dMatrix")
	n <- max(dim(Lr))
	Current <- .currentR(L, Lr, A, n, indexOrigin, indexGoal)
	result <- as(x,"RasterLayer")
	dataVector <- rep(0,times=ncell(result))
	dataVector[tc] <- Current
	result <- setValues(result, dataVector)
	return(result)
}

# Author: Jacob van Etten jacobvanetten@yahoo.com, based on Matlab code by Marco Saerens
# IE University
# Date :  January 2010
# Version 1.0
# Licence GPL v3

setMethod("passage", signature(x = "TransitionLayer", origin = "Coords", goal = "Coords", theta="numeric"), def = function(x, origin, goal, theta, totalNet="net", output="RasterLayer")
	{
		cellnri <- cellFromXY(x, origin)
		cellnrj <- cellFromXY(x, goal)
		x <- .transitionSolidify(x)
		tc <- transitionCells(x)

		ci <- match(cellnri,tc)
		cj <- match(cellnrj,tc)
		
		result <- .randomShPaths(x, ci, cj, theta, tc, totalNet, output)
		return(result)
	}
)

setMethod("passage", signature(x = "TransitionLayer", origin = "RasterLayer", goal = "RasterLayer", theta="numeric"), def = function(x, origin, goal, theta, totalNet="net", output="RasterLayer")
	{
		#check if Transition and RasterLayers coincide
		ci <- which(getValues(origin))
		cj <- which(getValues(goal))
		x <- .transitionSolidify(x)
		tc <- transitionCells(x)
		result <- .randomShPaths(x, ci, cj, theta, tc, totalNet, output)
		return(result)
	}
)

.randomShPaths <- function(x, ci, cj, theta, tc, totalNet, output)
{
	if(theta < 0 | theta > 20 ) {stop("theta value out of range (between 0 and 20)")}
	
	tr <- transitionMatrix(x,inflate=FALSE)
	
	trR <- tr
	trR@x <- 1 / trR@x 
	nr <- dim(tr)[1] 
	Id <- Diagonal(nr) 
	rs <- rowSums(tr)
	rs[rs>0] <- 1/rs[rs>0]
	P <- tr * rs

	W <- trR
	W@x <- exp(-theta * trR@x) #zero values are not relevant because of next step exp(-theta * trR@x) 
	W <- W * P 

	return(.probPass(x, Id, W, nr, ci, cj, tc, totalNet, output))
}
	
.probPass <- function(x, Id, W, nr, ci, cj, tc, totalNet, output)
{
	nc <- ncell(x)
	
	Ij <- Diagonal(nr)
	Ij[cbind(cj,cj)] <- 1 - 1 / length(cj)
	Wj <- Ij %*% W
	
	ei <- rep(0,times=nr)
	ei[ci] <- 1 / length(ci)
	
	ej <- rep(0,times=nr)

	ej[cj] <- 1 / length(cj)
	
	IdMinusWj <- as((Id - Wj), "dgCMatrix")
	
	zci <- solve(t(IdMinusWj),ei)
	zcj <- solve(IdMinusWj, ej)
	zcij <- sum(ei*zcj)
	
	if(zcij < 1e-300)
	{
		if(output == "RasterLayer")	
		{
			result <- as(x,"RasterLayer")
			result[] <- rep(0,times=nc)
		}	
		if(output == "TransitionLayer")
		{
			result <- x
			transitionMatrix(result) <- Matrix(0, nc, nc)
			result@transitionCells <- 1:nc
		}
	}
	
	else
	{
		# Computation of the cost dij between node i and node j
		# dij <- (t(zci) %*% (trR * Wj) %*% zcj) / zcij
	
		# Computation of the matrix N, containing the number of passages through
		# each arc
		N <- (Diagonal(nr, as.vector(zci)) %*% Wj %*% Diagonal(nr, as.vector(zcj))) / zcij
	
		#N is here the NET number of passages, like McRae-random walk

		if(output == "RasterLayer")
		{
			if(totalNet == "total")
			{
				# Computation of the vector n, containing the number of visits in
				# each node
				n <- pmax(rowSums(N),colSums(N)) #not efficient but effective
			}
			# Computation of the matrix Pr, containing the transition
			# probabilities
			#rn <- rep(0, times=length(n))
			#rn[n>0] <- 1 / n[n>0]
			#Pr <- N * rn
			#Pr <- N * (1 / n)
	
			#net visits
			if(totalNet == "net")
			{
				nNet <- abs(skewpart(N))
				n <- pmax(rowSums(nNet),colSums(nNet))
				n[c(ci,cj)] <- 2 * n[c(ci,cj)]
			}
			result <- as(x,"RasterLayer")
			dataVector <- rep(NA,times=nc)
			dataVector[tc] <- n
			result <- setValues(result, dataVector)
		}
	
		if(output == "TransitionLayer")
		{
			result <- x
			if(totalNet == "total")
			{
				tr <- Matrix(0, nc, nc)
				tr[tc,tc] <- N
			}
			if(totalNet == "net")
			{
				nNet <- skewpart(N) * 2
				nNet@x[nNet@x<0] <- 0
				tr <- Matrix(0, nc, nc)
				tr[tc,tc] <- nNet
			}
			
			transitionMatrix(result) <- tr
			result@transitionCells <- 1:nc
			
		}
	}
	return(result)
}