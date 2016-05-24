read.nex <- function(x){
	
	x <- scan(x, what = "c", quiet = TRUE)
		
	## eliminate comments
	## ------------------
	left <- grep("\\[", x)
	right <- grep("\\]", x)
	if ( length(left) > 0 ){
	  m <- cbind(left, right)
	  x <- x[-unlist(apply(m, 1, function(x) x[1]:x[2]))]
	}
	
	x <- x[x != ""]
	
  ## getting number of taxa
  ## ----------------------
	ntax <- x[grep("ntax", x, ignore.case = TRUE)]
	ntax <- gsub("[[:alpha:]]|[[:punct:]]", "", ntax )
	nb <- ntax <- as.numeric(unique(ntax))
		
	## getting number of characters	
  ## ----------------------------
	ncha <- x[grep("nchar", x, ignore.case = TRUE)]
	ncha <- gsub("[[:alpha:]]|[[:punct:]]", "", ncha )
	ncha <- as.numeric(unique(ncha))
	
	## get beginning and end of matrix
  ## -------------------------------
	start <- grep("^\t?matrix$", x, ignore.case = TRUE)
	end <- grep(";", x)
	end <- min(end[end > start])
	M <- x[(start + 1):(end - 1)]
	
	# assemble DNAbin object:
	# -----------------------
	nblock <- ceiling(ncha / nchar(M[2]))
	id <- seq(1, 2 * ntax, by = 2)
	nam <- M[id]
	fuse <- function(s, M, nblock, ntax){
	  paste(M[seq(s, length.out = nblock, by = ntax * 2)], collapse = "")
	}
	seq <- lapply(id + 1, fuse, M = M, nblock = nblock, ntax = ntax)
	obj <- list(nb = ntax, seq = seq, nam = nam, com = NA)
	class(obj) <- "alignment"
	as.DNAbin(obj)
}

