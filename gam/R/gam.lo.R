"gam.lo" <-
function(x, y, w = rep(1, length(y)), span = 0.5, degree = 1, ncols = p, xeval
	 = x)
{
	storage.mode(x) <- storage.mode(y) <- storage.mode(w) <- storage.mode(
		span) <- "double"
	storage.mode(degree) <- "integer"
	if(is.null(np <- dim(x))) {
		n <- as.integer(length(x))
		p <- as.integer(1)
	}
	else {
		np <- as.integer(np)
		n <- np[1]
		p <- np[2]
	}
	storage.mode(ncols) <- "integer"
	o <- gam.match(x)
	nef <- o$nef
	nvmax <- as.integer(200 + 300 * (1 - 1/log(max(c(nef - 200, 3)))))
	liv <- as.integer(50 + (2^ncols + 4) * nvmax + 2 * nef)
	lv <- as.integer(50 + (3 * ncols + 3) * nvmax + nef + (ifelse(degree ==
		2, ((ncols + 2) * (ncols + 1))/2, ncols + 1) + 2) * (nef * span +
		1))
	fit <- .Fortran("lo0",
		x,
		y,
		w,
		n,
		ncols,
		p,
		nvmax,
		span,
		degree,
		o$o,
		nef,
		df = double(1),
		s = double(n),
		var = double(n),
		beta = double(p + 1),
		iv = integer(liv),
		liv,
		lv,
		v = double(lv),
                integer(2*ncols),
		double(nef * (p + ncols + 8) + 2 * p + n + 9),
                        PACKAGE="gam")
	if(!missing(xeval)) {
		storage.mode(xeval) <- "double"
		m <- as.integer(dim(xeval)[1])
		if(length(m) == 0)
			m <- as.integer(length(xeval))
		.Fortran("lowese",
			fit$iv,
			liv,
			lv,
			fit$v,
			m,
			xeval,
			s = double(m),
                         PACKAGE="gam")$s - cbind(1, xeval) %*% fit$beta
	}
	else list(residuals = y - fit$s, var = fit$var, nl.df = fit$df)
}
