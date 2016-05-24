## PACKAGE: ips
## CALLED BY: USER, mrbayes
## AUTHOR: Christoph Heibl (at gmx.net)
## LAST UPDATE: 2014-08-11
write.nex <-
function(x, file = "", interleave = 60, taxblock = FALSE){
	
	if (class(x) == "alignment")
		stop ("Sequence alignment must be of class 'DNAbin'")
	
	str2cha <- function(x){unlist(strsplit(x, ""))}
	
	if (class(x) == "DNAbin") datatype <- "dna"
	if (class(x) == "dist") datatype <- "distances"
	if (class(x) == "data.frame") datatype <- "standard"
	missing <- if("?" %in% x[[1]]) "?" else "N"
	
	ntax <-  if (datatype == "distances") 
		attr(x, "Size") else dim(x)[[1]]
	nchar <- dim(x)[[2]]
	taxnames <- if (datatype == "distances") 
		labels(x) else rownames(x) 
	
	header <- c("#NEXUS", 
		paste("\n[created by ips on ", date(), "]\n", sep = ""))
	
	# TAXA BLOCK
	# -----------------
	if (taxblock){
		tb <- c("begin taxa;", 
			paste("\tdimensions ntax=", ntax,";", sep = ""),
			"\ttaxlabels",
			paste("\t", taxnames, sep = ""), ";\n")
	}
	
	# taxonnames of same length
	# -------------------------
	len <- nchar(taxnames)
	mlen <- max(len)
	len <- mlen - len + 1
	foo <- function(x){
		x <- rep(" ", x)
		paste(x, collapse = "")
	}
	ws <- lapply(len, foo)
	taxnames <- paste(taxnames, ws, sep = "")
	
	# indices of partitions
	# ---------------------
	int <- if (!is.numeric(interleave)) "" else " interleave"
	if (!interleave) interleave <- nchar
	nbpart <- ceiling(nchar/interleave)
	pt <- matrix(nrow = nbpart, ncol = 2)
	pt[1, ] <- c(1, interleave)
	if (nbpart > 1)
		for (i in 2:(dim(pt)[1])){
		    pt[i, ] <- c(pt[i - 1, 2] + 1, pt[i - 1, 2] + interleave)
		    pt[nbpart, 2] <- nchar
		}
		
	# assemble matrix
	# ---------------
	m <- "matrix"
	for (i in seq(along = pt[, 1])){
		sm <- as.character(x[, pt[i, 1]:pt[i, 2]])
		if (is.null(dim(sm))) sm <- as.matrix(sm, ncol = 1)
		sm <- apply(sm, 1, paste, collapse = "")
		sm <- paste(taxnames, sm)
		m <- c(m, sm)
	}
	m <- c(m, ";\nend;")
	
	
	# write DATA BLOCK
	# -----------------
	if (datatype == "dna" || datatype == "standard"){
	
		if (datatype == "standard" && !identical(x[1, 1], 
		    round(x[1, 1]))) 
			datatype <- "continuous"
		
		if (!datatype == "standard") dt <- datatype		
		else dt <- paste(datatype, " symbols=\"", 
		paste(unique(unlist(x)), collapse = " "), "\"", 		
		    sep = "")
		
	data <- if (taxblock) "characters" else "data"
	
	# assemble DATA BLOCK
	# -------------------
	db <- c(paste("begin ", data, ";", sep = ""),
		paste("\tdimensions ntax=", ntax, " nchar=", nchar, ";",
		    sep = ""),
		paste("\tformat datatype=", dt, " missing=", missing,
		    " gap=-", int, ";", sep = ""),
		m
	)
	}
	
	# assemble NEXUS file:
	# --------------------
	if (taxblock) nex <- c(header, tb, db)					
    else nex <- c(header, db)
		
	# write NEXUS file
	# ----------------
	if ( file == "" ) {
    cat(nex, sep = "\n")
    invisible(nex)
	} else {
    write(nex, file = file)
	}
}
