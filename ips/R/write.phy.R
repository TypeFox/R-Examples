write.phy <- function(x, file = "", interleave = FALSE, strict = FALSE){
	
	str2cha <- function(x){unlist(strsplit(x, ""))}
		
	datatype <- ifelse(is.numeric(x[1,1]), "continuous", "nc")
	
	ntax <- nrow(x)
	nchar <- ncol(x)
														
	# taxon names
	# -----------
	taxnames <- rownames(x)
	
	if ( strict ){
		taxnames <- substring(taxnames, 1, truncate)
		missing <- 10 - unlist(lapply(strsplit(taxnames, ""), length))
		for (i in seq(along = taxnames))
			taxnames[i] <- paste(taxnames[i], paste(rep("*",
			    missing[i]), collapse = ""), sep = "")
		if (any(duplicated(taxnames))) 
		    cat("WARNING: Truncation of taxon names created", 
		        "identical strings.")
	}
	else {
		xx <- nchar(taxnames)
		diff <- max(xx) - xx + 3 # add 3 for baseml
		for (i in 1:ntax) taxnames[i] <- paste(taxnames[i], 
			paste(rep(" ", diff[i]), collapse = ""), sep = "")
	}
	
	# indices of partitions: pt
	# -------------------------
	if (!interleave) interleave <- nchar
	nbpart <- ceiling(nchar/interleave)
	pt <- matrix(nrow = nbpart, ncol = 2)
	pt[1, ] <- c(1, interleave)
	if (nbpart > 1)
		for (i in 2:(dim(pt)[1])){
		    pt[i, ] <- c(pt[i - 1, 2] + 1, pt[i - 1, 2] + interleave)
		    pt[nbpart, 2] <- nchar
		}
	
	# assemble matrix: m
	# ------------------
	phy <- paste(ntax, nchar)
	for (i in seq(along = pt[, 1])){
		sm <- as.character(x[, pt[i, 1]:pt[i, 2]])
		if (is.null(dim(sm))) sm <- as.matrix(sm, ncol = 1)
		sm <- apply(sm, 1, paste, collapse = "")
		if (i == 1) sm <- paste(taxnames, sm)
		if (i < max(seq(along = pt[, 1]))) sm <- c(sm, "")
		phy <- c(phy, sm)
	}
	
	# write PHYLIP file
	# -----------------
	if (file == "") {
		cat(phy, sep = "\n")
	}
	else {
    write(phy, file = file)
	}
}
