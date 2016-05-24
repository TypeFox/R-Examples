combine.buffers <- function(...) {
	
	buf <- list(...)
	
	if (length(buf) == 0) {
		stop("You must provide at least one 'buffer' object.\n")
	} else { 
		if (length(buf) == 1) {
			if (!is(buf[[1]],'buffer')) stop("combine.buffers needs at least one 'buffer' object as produced by create.buffer().\n")
			warning("You provided only one buffer object: ther's nothing to combine.\n You will only get a new bathy object...\n")
			return(buf[[1]])
		} else {
			
			if (!all(sapply(buf, function(x) is(x,"buffer")))) stop("At least one object is not of class 'buffer'. See create.buffer() to produce such objects.")
			if (any(apply(sapply(buf, function(x) dim(x[[1]])), 1, function(x) length(table(x))) != 1)) stop("All buffer objects provided must contain a bathy object of compatible caracteristics (identical dimensions and location)")

			buf <- lapply(buf, function(x) x[[1]])

			lon <- rownames(buf[[1]])
			lat <- colnames(buf[[1]])

			# Stack buffers in a new array
			buf <- array(as.vector(sapply(buf,as.vector)),dim=c(dim(buf[[1]]),length(buf)))

			# Apply element wise function : if all elements in a position are NA, leave it that way, otherwise leave the depth/altitude value
			res <- apply(buf, 1:2, function(x) ifelse(any(!is.na(x)),na.omit(x)[1],NA))
			rownames(res) <- lon
			colnames(res) <- lat

			class(res) <- "bathy"

			return(res)
		}
	
	}
	
}