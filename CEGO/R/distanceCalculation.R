#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Calculate Distance Matrix
#'
#' Calculate the distance between all samples in a list, and return as matrix.
#'
#' @param X list of samples, where each list element is a suitable input for \code{distFun}
#' @param distFun Distance function of type f(x,y)=r, where r is a scalar and x and y are elements whose distance is evaluated.
#'
#' @return The distance matrix
#'
#' @examples
#' x <- list(5:1,c(2,4,5,1,3),c(5,4,3,1,2), sample(5))
#' distanceMatrix(x,distancePermutationHamming)
#'
#' @export
###################################################################################
distanceMatrix <-function(X,distFun){
	n <- length(X)
	m <- matrix(0,nrow=n, ncol=n)
	for(i in seq_len(n - 1))
		m[seq(i+1, n),i] <- m[i,seq(i+1, n)] <- distanceVector(X[[i]],X[seq(i+1, n)],distFun)
	m
}

###################################################################################
#' Calculate Distance Vector
#'
#' Calculate the distance between a single sample and all samples in a list.
#'
#' @param a A single sample which is a suitable input for \code{distFun}
#' @param X list of samples, where each list element is a suitable input for \code{distFun}
#' @param distFun Distance function of type f(x,y)=r, where r is a scalar and x and y are elements whose distance is evaluated.
#'
#' @return A numerical vector of distances
#'
#' @examples
#' x <- 1:5
#' y <- list(5:1,c(2,4,5,1,3),c(5,4,3,1,2))
#' distanceVector(x,y,distancePermutationHamming)
#'
#' @export
###################################################################################
distanceVector <-function(a,X,distFun){
	unlist(lapply(X,distFun,a))
}

###################################################################################
#' Update distance matrix
#'
#' Update an existing distance matrix \code{D_mat} by adding distances
#' of all previous candidate solutions to one new candidate solution, \code{d_vec= d(x_i,x_new)}.
#'
#' @param distanceMat original distance matrix \code{D_mat}
#' @param x list of candidate solutions, last in list is the new solution
#' @param distanceFunction Distance function of type f(x,y)=r, where r is a scalar and x and y are candidate solutions whose distance is evaluated.
#'
#' @return matrix of distances between all solutions x
#'
#' @examples
#' x <- list(5:1,c(2,4,5,1,3),c(5,4,3,1,2))
#' dm <- distanceMatrix(x,distancePermutationHamming)
#' x <- append(x,list(1:5))
#' dmUp <- distanceMatrixUpdate(dm,x,distancePermutationHamming)
#'
#' @export
#' @keywords internal
###################################################################################
distanceMatrixUpdate <- function(distanceMat,x,distanceFunction){
	count <- length(x)
	if(length(distanceFunction)==1){ # in case of a single distance function (all models)
		newdist = distanceVector(x[[count]],x[-count],distanceFunction)
		distanceMat = cbind(rbind(distanceMat,c(newdist)),c(newdist,0))
	}else{	# in case of multiple distance functions (kriging only atm.)
		for(i in 1:length(distanceFunction)){
			newdist = distanceVector(x[[count]],x[-count],distanceFunction[[i]])
			distanceMat[[i]] <- cbind(rbind(distanceMat[[i]],c(newdist)),c(newdist,0))
		}
	}
}

###################################################################################
#' Distance Matrix Wrapper
#'
#' Wrapper to calculate the distance matrix, with one or multiple distance functions.
#'
#' @param x list of candidate solutions whose distance is evaluated
#' @param distanceFunction Distance function of type f(x,y)=r, where r is a scalar and x and y are candidate solutions whose distance is evaluated.
#'
#' @return matrix of distances between all solutions in list x
#'
#' @examples
#' x <- list(5:1,c(2,4,5,1,3),c(5,4,3,1,2))
#' dm1 <- distanceMatrix(x,distancePermutationHamming)
#' dm2 <- distanceMatrix(x,distancePermutationInsert)
#' dmBoth <- distanceMatrixWrapper(x,list(distancePermutationHamming,distancePermutationInsert))
#'
#' @export
#' @keywords internal
###################################################################################
distanceMatrixWrapper <- function(x,distanceFunction){
	if(length(distanceFunction)==1){ # in case of a single distance function (all models)
		distances <- distanceMatrix(x,distanceFunction)
	}else{	# in case of multiple distance functions (kriging only atm.)
		distances <- list()
		for(i in 1:length(distanceFunction)){
			distances[[i]] <- distanceMatrix(x,distanceFunction[[i]]) 
		}
	}
  distances
}