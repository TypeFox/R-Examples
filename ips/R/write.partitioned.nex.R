write.partitioned.nex <-
function(..., file, labels){
	
  x <- list(...)
	sm <- cbind(..., fill.with.gaps = TRUE)
	sm <- sm[order(rownames(sm)),]
	ntax <- nrow(sm)
	nchar <- ncol(sm)
	
	# taxonnames of same length
  # -------------------------
	taxnames <- rownames(sm)
	len <- nchar(taxnames)
	mlen <- max(len)
	len <- mlen - len + 1
	fill.with.space <- function(x){
	  paste(rep(" ", x), collapse = "")
	}
	ws <- lapply(len, fill.with.space)
	taxnames <- paste(taxnames, ws, sep = "")
	
	
	# nucleotides in each partition:
	pt <- t(sapply(x, dim))
	pt[, 1] <- 1
	pl <- pt[, 2]
	for (i in 2:(dim(pt)[1])){
	  pt[i, ] <- pt[i, ] + pt[i - 1, 2]
	}
	
  ## header
  ## ------
	h <- c(
	  "#nexus",
	  paste("[created by ips on ", date(), "]\n", sep = ""),
	  "begin data;",
	  paste("\tdimensions ntax=", ntax, " nchar=", nchar, ";", sep = ""),
	  "\tformat datatype=dna missing=n gap=- interleave;"
    )

	## assemble matrix
	## ---------------
	m <- "matrix"
	
	for (i in seq(along = pt[, 1])){
		l <- paste("Position ", pt[i, 1], "-", pt[i, 2], ": ",  			
               labels[i], " (", pl[i], "bp)", sep = "")
		m <- c(m, paste("\n[", l, "]\n"))
		seq <- as.character(sm[, pt[i, 1]:pt[i, 2]])
		seq <- apply(seq, 1, paste, collapse = "")
		seq <- paste(taxnames, seq)
		m <- c(m, seq)
	}
	m <- c(m, ";\nend;")
	nex <- c(h, m)
	
	# write NEXUS file:
	# -----------------
	if ( missing(file) ) {
    nex
	}
	else {
    write(nex, file)
	}
}

