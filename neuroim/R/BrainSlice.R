#' @include AllClass.R
{}
#' @include AllGeneric.R
{}


#' BrainSlice constructor
#' @param data data vector or matrix
#' @param space an instance of class BrainSpace
#' @param indices linear indices corresponding used if \code{data} is a 1D vector.
#' @export
#' @examples
#' bspace <- BrainSpace(c(64,64), spacing=c(1,1))
#' dat <- array(rnorm(64*64), c(64,64))
#' bslice <- BrainSlice(dat,bspace)
#' print(bslice) 
BrainSlice <- function(data, space, indices=NULL) {
	if (ndim(space) != 2) {
		stop("incorrect dimension for BrainSlice")
	}
	
	if (is.null(indices)) {
		if (length(dim(data)) != 2) {
		  stopifnot(length(data) != prod(dim(space)[1:2]))
			data <- matrix(data, dim(space)[1], dim(space)[2])
		} 
		new("BrainSlice", .Data=data, space=space)
	
	} else {
		mdat <- matrix(0, dim(space))
		mdat[indices] <- data
		new("BrainSlice", .Data=mdat, space=space)
	}
}


#' @export
#' @rdname gridToIndex-methods
setMethod(f="gridToIndex", signature=signature(x = "BrainSlice", coords="matrix"),
		def=function(x, coords) {            
			dx <- dim(x)
			nsize <- prod(dx)
			apply(coords, 1, function(vox) {
						(vox[2]-1)*dx[1] + vox[1]
					})
		})



#' @export
#' @rdname indexToGrid-methods
setMethod(f="indexToGrid", signature=signature(x = "BrainSlice", idx="index"),
		def=function(x, idx) {            
			.indexToGrid(idx, dim(x))
		})






