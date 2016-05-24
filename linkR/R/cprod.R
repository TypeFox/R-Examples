cprod <- function(u, v, h='right'){

	if(is.matrix(u)){		
		if(nrow(u) == 1 || ncol(u) == 1) u <- as.vector(u)
	}
	if(is.matrix(v)){
		if(nrow(v) == 1 || ncol(v) == 1) v <- as.vector(v)
	}
	if(is.vector(v) && is.null(dim(v))){
		if(length(u) != length(v)){
			cat("CPROD ERROR: VECTORS ARE NOT OF EQUAL DIMENSION\n")
			return
		}
		if(length(v) == 2){
			u[3] <- 0
			v[3] <- 0
		}
		r <- c(u[2]*v[3] - u[3]*v[2], u[3]*v[1] - u[1]*v[3], u[1]*v[2] - u[2]*v[1])
		r <- unlist(r)
		names(r) <- names(u)
		if(h == 'right') return(r)
		return(-r)
	}
	if(is.null(dim(v))){
		r <- u[1]*v[2]-u[2]*v[1]
		if(h == 'right') return(r)
		return(-r)
	}
	if(is.matrix(v)){
		u <- unlist(u)
		r <- matrix(c(u[2]*v[3, ] - u[3]*v[2, ], u[3]*v[1, ] - u[1]*v[3, ], u[1]*v[2, ] - u[2]*v[1, ]), nrow=3, ncol=ncol(v), byrow=T)
		if(h == 'right') return(r)
		return(-r)
	}
	if(is.vector(v) || dim(v)[1] == 1 || dim(v)[2] == 1){
		r <- c(u[2]*v[3] - u[3]*v[2], u[3]*v[1] - u[1]*v[3], u[1]*v[2] - u[2]*v[1])
		r <- unlist(r)
		names(r) <- names(u)
		if(h == 'right') return(unlist(r))
		return(-unlist(r))
	}
}