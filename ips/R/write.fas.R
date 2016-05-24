write.fas <- function(x, file = "", interleave = FALSE, 
                      truncate = FALSE, append = FALSE){
  
  # x MUST be a list of class DNAbin
  # --------------------------------
  if ( is.matrix(x) ){
    x <- as.list(x)
  }
  if ( is.data.frame(x) ){
    y  <- as.list(x[, 1])
    names(y) <- rownames(x)
    x <- y
  }
		
	# taxonnames
	# ----------
	taxnames <- names(x)
	if (truncate){
		taxnames <- substring(taxnames, 1, truncate)
		if (any(duplicated(taxnames)))
			warning(paste("truncation of taxon names created",
			    "identical strings"))
	}
	taxnames <- paste(">", taxnames, sep = "")
	
	# function to interleave sequences:
	# ---------------------------------
	create.interleaved <- function(x, bp){
		x <- unlist(strsplit(x, ""))
		ncha <- length(x)
		if (!is.numeric(bp)) bp <- ncha
		nbpar <- ceiling(ncha / bp)
		out <- vector(length = nbpar)
		for (i in seq(along = out)){
			idmax <- i * bp
			if (idmax > ncha) idmax <- ncha
			idmin <- i * bp - bp + 1
			out[i] <- paste(x[idmin:idmax], collapse = "")
		}
		out
	}
	
	# assemble FASTA file:
	# --------------------
  s <- as.character(x)
  s <- lapply(s, paste, collapse = "")
	s <- lapply(s, create.interleaved, bp = interleave)
	nbpar <- sapply(s, length)
	s <- unlist(s)
	fas <- vector(length = length(s) + length(taxnames))
	id <- 1
	for (i in seq(along = nbpar[-1])) 
	    id <- c(id, tail(id, 1) + nbpar[i] + 1)
	fas[id] <- taxnames
	fas[-id] <- unlist(s)
	
	# write FASTA file
	# ----------------
	if ( missing(file) ) {
    cat(fas, sep = "\n")
	} else {
    write(fas, file = file, append = append)
	}
}
