NULL
#'
#' @export 
#' @rdname as.blockmatrix
#' 

as.blockmatrix <- function(M=NULL,...) {

	
	UseMethod("as.blockmatrix",M)
	
}







NULL
#print

#' @export 
#'  
#' @rdname as.blockmatrix
#' @method as.blockmatrix default
#' @S3method as.blockmatrix default
#' @aliases as.blockmatrix 
#' @export
#' 



as.blockmatrix.default <- function(M,adjust_zero=TRUE,zero_element="0",...) {
	
	
	out <- blockmatrix(...)
	out <- as.blockmatrix(out,adjust_zero=adjust_zero,add_zero_matrix=FALSE,zero_element=zero_element) 
	return(out)
	
}

NULL







NULL
#

#' @export 
#'  
#' @rdname as.blockmatrix
#' @method as.blockmatrix blockmatrix
#' @S3method as.blockmatrix blockmatrix
#' @aliases as.blockmatrix 
#' @export
#' 



as.blockmatrix.blockmatrix <- function(M,adjust_zero=TRUE,add_zero_matrix=FALSE,zero_element="0",...) {
	
	out <- M 
	if (is.zero.blockmatrix(out)) return(out)
	
	if (adjust_zero) {
		
		nrow <- nrow(M)
		ncol <- ncol(M)
		value <- value(M)
		value0 <- value
	
		for (i in 1:nrow) {
			for (j in 1:ncol) {
				temp <- out[i,j]
				
				if (length(temp)==1) {
					if (temp==0) value[i,j]=zero_element
					
				} else if (length(temp)>1) {
					
					mx <- max(temp,na.rm=TRUE)
					mn <- min(temp,na.rm=TRUE)
			
					
					if ((mn==0) & (mx==0)) value[i,j]=zero_element
				}
			}
		}

		if (min(value==zero_element)==1) {
			
			out <- NULL
			out$value <- 0
			class(out) <- "blockmatrix"
			
		} else {
		
			names <- as.vector(value[value==value0])
			
            l <- M[names]		
			l <- l[!is.na(names(l))]
			
			out <- blockmatrix(list=l,value=value,use.as.blockmatrix=FALSE) ## Using use.block.matrix=TRUE returns Error: evaluation nested too deeply: infinite recursion / options(expressions=)?
		
		}
		
		
	} else if (add_zero_matrix) {
		
		nrowe <- nrow_elements(M,zero_element=zero_element)
		ncole <- ncol_elements(M,zero_element=zero_element)
		
		if ((!is.na(nrowe)) & (!is.na(ncole))) M[zero_element] <- array(0,c(nrowe,ncole)) 
		
	} else {
		out <- M 
	}

#	value <- value(out)
#	list <- out[as.vector(value(out))]
	
	
	
	return(out)
	
}

NULL

#' 
#'\code{as.blockmatrix} S3 method for \code{blockmatrix}, \code{matrix} and \code{NULL} object
#'
#' @param M a \code{matrix} or \code{blockmatrix} object
#' @param nrowe number of rows for each block (element of the blockmatrix)
#' @param ncole number of columns for each block (element of the blockmatrix)
#' @param nrow number of rows for block-matrix 
#' @param ncol number of columns of blockmatrix
#' @param adjust_zero logical value. If \code{TRUE} (Default) it replaces the zero matrices with \code{zero_element}.
#' @param add_zero_matrix logical value.  If \code{TRUE} it adds a zero-element element matrix as an object called \code{zero_element} in the blockmatrix
#' @param zero_element see \code{\link{ncol_elements}} or \code{\link{nrow_elements}}
#' @param ... further arguments 
#' @export
#' @rdname as.blockmatrix
#' @method as.blockmatrix matrix
#' @S3method as.blockmatrix matrix
#' @aliases as.blockmatrix
#'
#' @author Emanuele Cordano 
#' 
#' 
# CREARE FUNZIONE DA FARE GIOVEDI!!!! 
as.blockmatrix.matrix <- function (M,nrowe=2,ncole=2,nrow=NULL,ncol=NULL,adjust_zero=TRUE,zero_element="0",...)  {
	
	if (is.null(nrow)) nrow=nrow(M)/nrowe
		
	if (is.null(ncol)) ncol=ncol(M)/ncole

	nrowe <- nrow(M)/nrow
	ncole <- ncol(M)/ncol
	
	list <- list()
	
	value <- array("V",c(nrow,ncol))

	for (c in 1:ncol) {
		for (r in 1:nrow) {
			rows <- ((r-1)*nrowe+1):(r*nrowe)
			cols <- ((c-1)*ncole+1):(c*ncole)		
	
		
			value[r,c] <- paste(value[r,c],r,",",c,sep="")
			list[[value[r,c]]] <- M[rows,cols]
			
		}
	}
	
	out <- as.blockmatrix(value=value,list=list,adjust_zero=adjust_zero,zero_element=zero_element)
	
	
	
	
	return(out)
	
	
}