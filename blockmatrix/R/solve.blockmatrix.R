NULL


#' 
#' \code{dim} S3 solve for \code{blockmatrix} object
#' 
#' @param a a \code{blockmatrix} or numeric object 
#' @param b a \code{blockmatrix} or numeric object. If omitted, it is \code{NULL}. See Details. 
#' @param as.blockmatrix logical value. If \code{TRUE} (Default), the output is a \code{blockmatrix} object
#' @param ... further arguments for method \code{solve}
#' 
#' @export
#' @title solve
#' @rdname solve
#' @method solve blockmatrix
#' @S3method solve blockmatrix
#' @aliases solve
#'
#' @author Emanuele Cordano 
#' @return the object \code{x} such that \code{a * x = b} where \code{*} is the matrix product.
#' 
#' @note If \code{b} is missing, i.e. \code{NULL}, it will be replaced by the corresponding identity matrix. So \code{x} is calculated as the right inverse of \code{a}.
#' The matrix system must be nonsingular and nonhomogeneous.




solve.blockmatrix <- function (a,b=NULL,as.blockmatrix=TRUE,...)  {
	
	
	nrowe <- ncol_elements(a)
	am <- as.matrix(a)
	
	if (!is.null(b)) {
		ncole <- ncol_elements(b)
		bm <- as.matrix(b)
	} else {
		ncole <- nrow_elements(a)
		bm <- diag(nrow(am))
	}
	
	am <- as.matrix(a)
	
	oo <- solve(am,bm,...)
	if ((is.null(nrowe)) | (is.null(ncole)) | (is.na(nrowe)) | (is.na(ncole))) as.blockmatrix=FALSE 
	if (as.blockmatrix) {
		out <- as.blockmatrix(oo,nrowe=nrowe,ncole=ncole)
	} else {
		out <- oo
		
	}
	return(out)
	
	
	
	
	
	
	
}