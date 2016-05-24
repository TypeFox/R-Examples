# gives formula terms as a list of symbols (interactions as sub-list)
# [note that it does not expand formulas]
# termlist(terms(~a * b+ c, simplify = TRUE))
## termlist(~a+b+a:b) --> list(a, b, list(a, b))
termlist <- function(x) {
	is.plus <- function(x) is.call(x) && x[[1L]] == "+"
	## parses interaction expression into list: a:b:c --> list(a,b,c)
	intr <- function(x) {
		# is it an expression for interaction? (e.g. a:b:c)
		is.intr <- function(x) is.call(x) && x[[1L]] == ":"
		if(is.intr(x)) {
			res <- list()
			repeat {
				res <- c(x[[3L]], res)
				x <- x[[2L]]
				if(!is.intr(x)) break
			}
			list(c(x, res))
		} else x
	}
	if(x[[1L]] == "~") x <- x[[length(x)]]
	res <- list()
	while(is.plus(x)) {
		res <- c(intr(x[[3L]]), res)
		x <- x[[2L]]
	}
	res <- c(intr(x), res)
	res	
}

# calculates all lower order term names:
# expandintr(1:3) --> c("1", "2", "1:2", "3", "1:3", "2:3", "1:2:3")
expandintr <- function(x) {
	asstr <- function(x) asChar(x, backtick = TRUE)
	if(!is.language(x)) {
		a <- sapply(x, asstr)
		k <- length(a)
		vapply(seq.int(2L^k - 1L), function(y) paste0(a[as.logical(intToBits(y)[1L:k])],
			collapse = ":"), "")
	} else asstr(x)
}

# given a formula, 'term dependency matrix', i.e. dependency of higer
# order terms on other lower order terms
termdepmat <- function(f) {
	trm <- terms(f, simplify = TRUE)
	tl <- termlist(trm)
	v <- attr(trm, "term.labels")
	n <- length(v)
	mat <- matrix(FALSE, n, n, dimnames = list(v, v))
	for(i in seq.int(n)) mat[match(expandintr(tl[[i]]), v), i] <- TRUE
	mat
}

# alternative to 'termdepmat', gives matrix dimension names as numbers
# so a,b,a:b  become 1,2,1:2 
termdepmat2 <- function(f) {
	filist <- formula2idx(f, asCall = FALSE)
	n <- length(filist)
	v <- vapply(filist, paste0, "", collapse = ":")
	mat <- matrix(FALSE, n, n, dimnames = list(v, v))
	for(i in seq.int(n)) mat[match(expandintr(filist[[i]]), v), i] <- TRUE
	mat
}

## combines term-dependency-matrices
#termdepmat_list <- function(fl) 
#	termdepmat_combine(lapply(fl, termdepmat))

termdepmat_combine <- function(x) {
	dm <- sum(vapply(x, nrow, 1L))
	mat <- matrix(FALSE, dm, dm)
	j <- 1L
	for(i in seq_along(x)) {
		n <- nrow(x[[i]])
		k <- seq.int(j, length.out = n)
		mat[k, k] <- x[[i]]
		j <- j + n
	}
	dn <- unlist(lapply(x, rownames))
	dimnames(mat) <- list(dn, dn)
	mat
}

# converts formula to a(n unevaluated) list of numeric indices
# e.g. a*b --> list(1,2,1:2)
formula2idx <- function(x, asCall = TRUE) {
	if(!is.call(x) || !inherits(x, "formula")) stop("'x' is not a formula")
	fac <- attr(delete.response(terms(x)), "factors")
	dimnames(fac) <- NULL
	ret <- apply(fac > 0L, 2L, which)
	if(asCall) as.call(c(as.name("list"), ret)) else ret 
}

formula_margin_check <- function(j, m) {
	stopifnot(is.logical(j))
	!any(m[!j, j], na.rm = TRUE)
}