print.xyz <- function(x, ...) {

	## Print a summary of bio3d 'xyz' object features
  	if(!inherits(x, "xyz")) {
    	stop("Input should be a bio3d 'xyz' object")
  	}

  	if( is.null(nrow(x)) )
  		x <- t(as.matrix(x))

	cat( paste0("\n   Total Frames#: ", nrow(x), 
		        "\n   Total XYZs#:   ", ncol(x), 
		        ",  (Atoms#:  ",  round(ncol(x)/3,3), 
		        ")\n\n") ) 

	if(ncol(x) > 7) {
	s <- paste("    [1]  ", paste(round(x[1,1:3],3),collapse="  "), "  <...>  ",
		       paste(  round(x[nrow(x),(ncol(x)-2):ncol(x)], 3), collapse="  "), 
		       "  [",length(x),"]", sep="")
	} else {
		s <- paste("    [1]  ", paste(round(x[1,],3),collapse="  "), 
				   "  [",length(x),"]", sep="")
	}
	cat(s,"\n\n")

	i <- paste(attributes(x)$names, collapse = ", ")
	j <- paste("Matrix DIM =", nrow(x), "x", ncol(x))
	cat(strwrap(paste(" + attr:", i, "\n",j), width = 45, exdent = 8),
        sep = "\n")
}

