##
##  s t r f i n d . R  Find Substrings
##


strfind <- function(s1, s2, overlap = TRUE) {
	stopifnot(is.character(s1), is.character(s2))
	if (length(s2) != 1)
		stop("Pattern 's2' must be a character vector of length 1.")

	L <- list()
	for (i in 1:length(s1)) {
		si <- s1[i]
		if (nchar(si) < nchar(s2)) {
			L[[i]] <- c()
		} else {
			L[[i]] <- findstr(s2, si, overlap=overlap)
		}
	}
	if (length(s1) == 1) L <- unlist(L)
	return(L)
}

strfindi <- function(s1, s2, overlap = TRUE) {
    stopifnot(is.character(s1), is.character(s2))

    strfind(tolower(s1), tolower(s2), overlap=overlap)
}

findstr <- function(s1, s2, overlap = TRUE) {
	stopifnot(is.character(s1), is.character(s2))
	if (length(s1) != 1 || length(s2) != 1)
		stop("Arguments must be character vectors of length 1.")
	if (nchar(s1) > nchar(s2)) {
		s <- s1; s1 <- s2; s2 <- s
	}
	if (s1 == '') return(c())

	n1 <- nchar(s1)
	n2 <- nchar(s2)
	r  <- c()
	i  <- 1
	while (i <= n2 - n1 + 1) {
		if (s1 == substring(s2, i, i  + n1 - 1)) {
			r <- c(r, i)
			if (!overlap) i <- i + n1 - 1
		}
		i <- i + 1
	}
	return(r)
}
