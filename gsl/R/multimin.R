# An R wrapper for GSL's multimin family of functions for minimizing functions.
#
# Written in 2007 by Andrew Clausen <clausen@econ.upenn.edu>
# - Added Nelder-Mead convergence interface in 2010.

#dyn.load("gsl.so")

multimin <- function(..., prec=0.0001)
{
	is.converged <- function(state, old.x)
	{
		if (state$method == "nm")
			return(multimin.fminimizer.size(state) < prec)

		convergence <- ifelse(state$is.fdf, state$df, old.x - state$x)
		return (sum(abs(convergence)) < prec)
	}

	state <- multimin.init(...)
	old.x <- state$x
	while (TRUE)
	{
		state <- multimin.iterate(state)
		if (is.converged(state, old.x))
			break
		old.x <- state$x
	}
	state
}

multimin.init <- function(x, f, df=NA, fdf=NA, method=NA, step.size=NA, tol=NA)
{

        multimin.method.names <-
          list(
               "conjugate-fr",     # fletcher-reeves
               "conjugate-pr",     # polak-ribiere
               "bfgs",             # broyden-fletcher-goldfarb-shanno
               "steepest-descent", 
               "nm"                # nelder-mead
               )

        multimin.method.f <- c(FALSE, FALSE, FALSE, FALSE, TRUE)
        multimin.method.fdf <- c(TRUE, TRUE, TRUE, TRUE, FALSE)

	stopifnot(length(formals(f)) == 1)
	stopifnot(is.numeric(x))

	n <- length(x)

	is.fdf = !missing(df)
	if (is.fdf) {
		stopifnot(length(formals(df)) == 1)
		if (is.na(fdf))
			fdf <- function(x) list(f=f(x), df=df(x))
		stopifnot(length(formals(fdf)) == 1)
		fdf_ <- function(x) {
			result <- new.env()
			vals <- fdf(x)
			result$f <- vals$f
			result$df <- vals$df
			result
		}

		if (missing(method))
			method <- "bfgs"
		if (missing(step.size))
			step.size <- 1
		if (missing(tol))
			tol <- 1
		stopifnot(is.numeric(tol))
	} else {
		if (missing(step.size))
			step.size <- rep(1, n)
		if (missing(method))
			method <- "nm"
	}

	stopifnot(is.numeric(step.size))

	if (is.character(method)) {
		method.name <- method
		method <- match(method, multimin.method.names)
		if (is.na(method))
			stop(paste(c(
		"The optimization method '", method.name,
		"' is not an option.  Try one of these:\n",
		paste(multimin.method.names, collapse=", "),
		".",
		sep="")))
	}
	method <- as.integer(method)
	stopifnot(1 <= method && method <= 5)
	if (!is.fdf && !multimin.method.f[[method]])
		stop(paste(
		"The optimization method '", multimin.method.names[[method]],
		"' needs a derivative function.\nIf you don't want to provide ",
		"one, these methods don't need derivatives:\n",
		paste(subset(multimin.method.names, multimin.method.f),
		      collapse=", "),
		".",
		sep=""))
	if (is.fdf && !multimin.method.fdf[[method]])
		stop(paste(
		"The optimization method '", multimin.method.names[[method]],
		"' doesn't use derivatives.\nIf you want to exploit a ",
		"derivative, use one of these methods:\n",
		paste(subset(multimin.method.names, multimin.method.fdf),
		      collapse=", "),
		".",
		sep=""))

	internal.state <- new.env()
	internal.state$f <- body(function(x) f(x))
	internal.state$n <- n
	internal.state$rho <- new.env()
	if (is.fdf) {
		internal.state$df <- body(function(x) df(x))
		internal.state$fdf <- body(fdf_)
	}

	if (is.fdf) {
		.Call("multimin_fdf_new", internal.state, x, method, step.size,
		      tol)
	} else {
		.Call("multimin_f_new", internal.state, x, method, step.size)
	}

	list(internal.state = internal.state, x=x, f=NA, df=rep(NA, n),
	     is.fdf=is.fdf, method=multimin.method.names[[method]])
}

multimin.iterate <- function(state)
{
	internal.state <- state$internal.state
	if (state$is.fdf) {
		state$code <- .Call("multimin_fdf_iterate", internal.state)
		state$x <- .Call("multimin_fdf_state_argmin", internal.state)
		state$f <- .Call("multimin_fdf_state_min", internal.state)
		state$df <- .Call("multimin_fdf_state_grad", internal.state)
	} else {
		state$code <- .Call("multimin_f_iterate", internal.state)
		state$x <- .Call("multimin_f_state_argmin", internal.state)
		state$f <- .Call("multimin_f_state_min", internal.state)
	}
	state
}

multimin.restart <- function(state)
{
	if (state$is.fdf)
		.Call("multimin_restart", state$internal.state)
	state
}

# Convergence criterion for Nelder-Mead
multimin.fminimizer.size <- function(state)
{
	stopifnot(!state$is.fdf)
	.Call("multimin_fminimizer_size", state$internal.state)
}

