B <- function(..., knots = NULL, knot.range = 0:5) {
    vars <- as.list(substitute(list(...)))[-1]
	if(length(vars) == 0) 
		stop("B() is not supported")
    if(length(vars) > 1)
        stop("B() with more than one variable is not supported")
	if(is.null(knots) && is.null(knot.range))
		stop("either 'knots' or 'knot.range' must be specified with B()")
    
    knots <- as.vector(knots)
    if(!is.null(knots)) {
        if(!is.numeric(knots))
            stop("knots must be numeric")
        if(length(knots) != length(unique(knots)))
            stop("knots must be unique")
        knots <- sort(knots)
        knot.range <- length(knots)
    }

	term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
	if (term[1] == ".")
		stop("B(.) is not supported")
	termname <- attr(terms(reformulate(term[1])),"term.labels")
	termlabel <- paste("B(",termname,")",sep="")
	ret <- list(term = term, termlabel = termlabel, knots = knots, 
                    knot.range = knot.range)
    class(ret) <- "B.smooth"
    ret
}