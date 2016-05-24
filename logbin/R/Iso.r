Iso <- function(...) {
    vars <- as.list(substitute(list(...)))[-1]
	if(length(vars) == 0) 
		stop("Iso() is not supported")
    if(length(vars) > 1)
        stop("Iso() with more than one variable is not supported")

	term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
	if (term[1] == ".")
		stop("Iso(.) is not supported")
	termname <- attr(terms(reformulate(term[1])),"term.labels")
	termlabel <- paste("Iso(",termname,")",sep="")
	ret <- list(term = term, termlabel = termlabel, knots = NA, 
                    knot.range = NA)
    class(ret) <- "Iso.smooth"
    ret
}