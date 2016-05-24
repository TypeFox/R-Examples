#' Creates the inverse covariance matrix for an intrinsic conditionally
#' autoregressive spatial model.
#' 
#' This function creates the ICAR precision matrices used in the spatial models
#' 
#' Constructs the inverse covariance matrix (aside from scaling) for the ICAR
#' model
#' 
#' @param xy An n x 2 matrix of spatial coordinates
#' @param threshold Distance threshold for neighborhood definition
#' @param rho The autoregressive parameter. Defaults to 1, which is the
#' Intrinsic Conditionally AutoRegressive model (ICAR)
#' @param fun If \code{TRUE} this function returns a function of rho that generates the precision matrix of a ICAR process
#' @return An n x n matrix
#' @author Devin S. Johnson <devin.johnson@@noaa.gov>
#' @export
#' @import Matrix
icar.Q <-
function(xy, threshold, rho=1, fun=FALSE)
{
	thresh <- threshold
	distmat <- as.matrix(dist(cbind(xy[,1], xy[,2])))
	dist.ind <- (distmat <= thresh)*1.0
	diag(dist.ind) <- 0
	num <- as.vector(apply(dist.ind,1,sum))
	D <- diag(num)
	if(fun) return(function(x){Matrix(D-x*dist.ind)})
	else return(D-rho*dist.ind)
}

