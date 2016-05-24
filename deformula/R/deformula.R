
deformula.zeroinf <- function(f, ..., zero.eps = 1.0e-12,
	rel.tol = 1.0e-8, start.divisions = 8, max.iter = 12) {
	ff <- function(x) f(x, ...)
	res <- .Call(integrate_zero_to_inf, body(ff), environment(),
		zero.eps, rel.tol, start.divisions, max.iter)
	names(res) <- c("value", "x", "w", "t", "h", "message")
    switch(as.character(res$message),
		"0"={res$message <- "OK"},
		"2"={res$message <- "Error: Some values become NaN."},
        stop("Unknown error code.")
    )
	res
}

deformula.moneone <- function(f, upper, lower, ..., zero.eps = 1.0e-12,
	rel.tol = 1.0e-8, start.divisions = 8, max.iter = 12) {
	ff <- function(x) {
		(upper - lower) * f(((upper - lower) * x + (upper + lower)) / 2.0, ...) / 2.0
	}
	res <- .Call(integrate_mone_to_one, body(ff), environment(),
		zero.eps, rel.tol, start.divisions, max.iter)
	names(res) <- c("value", "x", "w", "t", "h", "message")
	res$x <- ((upper - lower) * res$x + (upper + lower)) / 2.0
    switch(as.character(res$message),
		"0"={res$message <- "OK"},
		"2"={res$message <- "Error: Some values become NaN"},
        stop("Unknown error code.")
    )
	res
}
