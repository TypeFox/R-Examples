##############################################
##checks whether data are of required format##
##############################################
check.data <-
function(X, format = c("probability.different", "percent.same", "general")){
	# X: the data
	# format: the data format to be checked; the default is "probability.different"

	format <- match.arg(format)
	if ((data.class(X) != "matrix") && (data.class(X) != "data.frame"))
		return(message("data must be matrix or data frame"))
    if (dim(X)[1] != dim(X)[2])
    	return(message("data must have same number of rows and columns"))
	M <- as.matrix(X)
  	if (mode(M) != "numeric")
		return(message("data must be numbers"))
	if (any(!is.finite(M)))
		return(message("data must be finite, i.e., not be NA, NaN, Inf, or -Inf"))
	# preceding queries constitute the "general" data format
	# additional queries for the "probability-different" and "percent-same" data formats are:
	if (format == "probability.different"){
		if (any(M < 0) || any(M > 1))
			return(message("data must be greater-equal 0 and less-equal 1"))
	}
	if (format == "percent.same"){
		if (any(M < 0) || any(M > 100))
			return(message("data must be greater-equal 0 and less-equal 100"))
	}

	# labeling of the rows and columns
	n <- dim(M)[1]
	permute.letters <- character(length(letters))  # for permuting the elements of the character vector 'letters = (a, b, c, ..., z)'
	permute.letters[1] <- letters[length(letters)]  # first element is set 'z'
	for (i in 1:25)  # second element is set 'a', third 'b', etc.
		permute.letters[i + 1] <- letters[i]
	# (1) if 'dimnames == NULL' or both 'rownames==colnames==NULL', and (1.1) if there are at most 234 objects, 
	#     the object names are set a1, b1, ..., z1, a2, b2, ..., z2, ..., etc., up to, if necessary,
	#     a9, b9, ..., z9 (i.e., the matrix may be as large as 234-times-234)
	if ((is.null(dimnames(M))) || (is.null(rownames(M)) && is.null(colnames(M)))){
		if (n <= (length(letters) * 9)){
		rownames(M) <- character(n)
		colnames(M) <- character(n)
		for (i in 1:n)
			rownames(M)[i] <- colnames(M)[i] <- paste(permute.letters[(i %% 26) + 1], ceiling(i / 26), sep = "")
		} else # and (1.2) if there are more than 234 objects, the object names are set v1, v2, ..., vn
			rownames(M) <- colnames(M) <- paste("v", 1:n, sep = "")
	}
	# (2) if colnames exist and rownames do not; colnames are assigned as rownames
	if (is.null(rownames(M)) && !is.null(colnames(M)))  
		rownames(M) <- colnames(M)
	# (3) if rownames exist and colnames do not; rownames are assigned as colnames
	if (!is.null(rownames(M)) && is.null(colnames(M)))  
		colnames(M) <- rownames(M) 
	# (4) if both rownames and colnames exist but do not coincide; produces a respective message
	if (!identical(rownames(M), colnames(M)))
		return(message("if labeling of both the rows and columns is provided manually, the same labeling must be used"))

	# checking duplicate names
	for (i in 1:n){
		if (charmatch(colnames(M)[i], colnames(M)) == 0)
			return(message(paste("duplicate name is found: ", colnames(M)[i], sep = "")))
	}

	return(M)
}
