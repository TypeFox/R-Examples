NULL

#' 
#'\code{Math} and \code{Ops} group of S3 methods for \code{blockmatrix} object
#'
#' @param x,e1,e2  \code{blockmatrix} objects 
#' @param as.blockmatrix logical value. If \code{TRUE} (Default), the output is a \code{blockmatrix} object
#' @param ... further arguments 
#' 
#'  @export
#' @rdname Math
#' @method Math blockmatrix
#' @S3method Math blockmatrix
#' @aliases Math
#'
#' @author Emanuele Cordano 
#' 
#' 


Math.blockmatrix <- function (x,as.blockmatrix=TRUE,...)  {
	

	T <- as.matrix(x)
	ncole=ncol_elements(x)
	nrowe=nrow_elements(x)
	
	if ((is.na(nrowe) | (is.null(nrowe)) | (is.na(ncole)) | (is.null(ncole)))) as.blockmatrix=FALSE

	O <- eval(call(.Generic,T))
	if (as.blockmatrix) 
		out <- as.blockmatrix(O,ncole=ncol_elements(x),nrowe=nrow_elements(x))
	else {
		out <- O
	}
	
	return(out)
}

NULL
#' @export
#' @rdname Math
#' @method Ops blockmatrix
#' @S3method Ops blockmatrix
#' @aliases Ops


Ops.blockmatrix <- function (e1,e2) {
	
	out <- NULL
	
	if (is.zero.blockmatrix(e1)) e1 <- 0 
	if (is.zero.blockmatrix(e2)) e2 <- 0 
	
	
	ncole <- ncol_elements(e2)

	
	

	if (class(e1)=="blockmatrix") {
		if (is.zero.blockmatrix(e1)) {
			em1 <- 0
		} else {
			em1 <- as.matrix(e1)
		}
	} else {
		em1 <- e1
#		nrowe <- nrow_elements(M)
	}
	
	if (class(e2)=="blockmatrix") {
		if (is.zero.blockmatrix(e2)) {
			em2 <- 0
				
		} else {
			em2 <- as.matrix(e2)
		}
			
		
		} else {
			em2 <- e2
#			ncole <- ncol_elements(e2)
	}
	nrowe <- nrow_elements(e1)
	ncole <- ncol_elements(e2)
	if ((is.na(nrowe)) | (is.null(nrowe))) nrowe <- nrow_elements(e2)
	if ((is.na(ncole)) | (is.null(ncole))) ncole <- ncol_elements(e1)

	
	oo <- eval(call(.Generic,em1,em2))
	out <- as.blockmatrix(oo,nrowe=nrowe,ncole=ncole)
	return(out)
}
