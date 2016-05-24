# ============================================================
# Low level functions to make special kinds of sparse matrices
# ============================================================

# make type-token (tt) Matrix from vector
# result: Types x Tokens, types are returned as separate rownames

ttMatrix <- function(vector, collation.locale = "C", simplify = FALSE) {

  # just is case the complete vector is NA, return NULL
  # which might occur in higher-level functions working on large datasets
  if (sum(!is.na(vector)) == 0) {
    return(NULL)
  } else {
    
  	# change locale for collation, defaults to pure unicode locale "C"
  	# setting NULL takes current locale on system
  	Sys.getlocale("LC_COLLATE") -> current.locale
  	if (!is.null(collation.locale)) {
  		Sys.setlocale("LC_COLLATE", collation.locale)
  	}
  
  	# factorization
  	factor <- factor(vector) # remove non-used levels
  	indices <- as.numeric(factor)
  	names <- levels(factor)
  
  	# just in case that there is missing data
  	indices <- na.omit(indices)
  	available <- which(!is.na(factor))
  	rows <- max(indices)
  	cols <- length(factor)
  
  	# make sparse matrix
  	M <- sparseMatrix(	i = indices,
  						j = available,
  						dims = c(rows,cols)
  						)
  
  	# change locale back to original
  	Sys.setlocale("LC_COLLATE", current.locale)
  	
  	if (simplify) {
  		rownames(M) <- names
  		colnames(M) <- vector
  		return(M)
  	} else {			
  		return(	list(	M = M, 
  						rownames = names
  						))
  	}
  }
}

# make part-whole (pw) Matrix from tokenized strings
# result: Segments x Strings, segments are returns as separate rownames
# gap is needed for not obtaining overlap for ngrams.
# gap.length = 1 is sufficient for bigrams, gap.length = 2 for 3-grams, etc.

pwMatrix <- function(strings, sep = "", gap.length = 0, gap.symbol = "\u00B7", simplify = FALSE) {

	# just to be sure that we are dealing with strings
	strings <- as.character(strings)
	
	# split the strings by specified separator
	parsed.strings <- strsplit(strings,split=sep)
	
	# count the number of segments per string	
	strings.length <- sapply(parsed.strings,length)
	
	# the following is simply used when no gaps are necessary
	# also a special case when there is only one string, just to catch errors
	if (gap.length == 0 | length(strings) == 1) {
		segments <- unlist(parsed.strings)
		indices <- 1:length(segments)		
	}

	# for unconnected lists of strings, add gaps to get ngrams right later on
	if (gap.length > 0 & length(strings) > 1) {
		
		# add gaps after strings and make one long vector with all segments
		gap <- rep(gap.symbol, gap.length)
		add.gap <- function(x) {c(x,gap)}
		segments <- unlist(sapply(parsed.strings, add.gap, simplify = FALSE))
		
		# and remove gap at the end of the long string
		segments <- head(segments,-gap.length)
		
		# make indices, and ignore gaps for segment indices in matrix
		indices <- (1:length(segments))[-which(segments == gap.symbol)]		
	}
						
	# part-whole Matrix: segments x strings
	M <- sparseMatrix(	i = indices,
						j = rep.int(1:length(strings),strings.length),
						dims = c(length(segments), length(strings))
						)

	if (simplify) {
		rownames(M) <- segments
		colnames(M) <- strings
		return(M)
	} else {			
		return(	list(	M = M, 
						rownames = segments
						))
	}
}


# Harmonize (alike to SQL "join") two matrices on a dimensions that have the same entities, 
# but in a different order (and possibly with different subsets)
# The idea is to take two factors, and return two Type-Token matrices, 
# in which the types (in the rows) are harmonized, and returned as rownames
# The matrix t(M1) %*% M2 (with the harmonized rows in the middle) can be used to JOIN two tables.

jMatrix <- function(rownamesX, rownamesY, collation.locale = "C") {

	# joined matrix
	J <- ttMatrix(c(rownamesX,rownamesY), collation.locale =  collation.locale)
	rownames <- J$rownames

	# split the joined matrix
	M <- t(J$M)
	M1 <- t(head(M,length(rownamesX)))
	M2 <- t(tail(M,length(rownamesY)))
	
	return(	list(	M1 = M1,
					M2 = M2,
					rownames = rownames
					))
}

jcrossprod <- function(X, Y, rownamesX = rownames(X), rownamesY = rownames(Y)) {
	
	J <- jMatrix(rownamesX, rownamesY)

	if (is(X,"nMatrix") & is(Y,"nMatrix")) {
		M <- crossprod( J$M1 %*% X, J$M2 %*% Y )
	} else {
		M <- crossprod( (J$M1*1) %*% X, (J$M2*1) %*% Y )
	}
	return(M)
}

tjcrossprod <- function(X, Y, colnamesX = colnames(X), colnamesY = colnames(Y)) {
	
	J <- jMatrix(colnamesX, colnamesY)

	if (is(X,"nMatrix") & is(Y,"nMatrix")) {
		M <- tcrossprod( X %*% t(J$M1), Y %*% t(J$M2) )
	} else {
		M <- tcrossprod( X %*% t(J$M1*1), Y %*% t(J$M2*1) )
	}
	return(M)
}
