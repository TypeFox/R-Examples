	plotclean <- function(x, y, id=NULL, data=parent.frame(), n=length(x), par.out=list(pch=20), ...) {
#	plot growth curves to identify outlying points and curves
#	plot y ~ x by id with data
#	identify up to n points
#	par.out passes graphical par args for outliers
	mcall <- match.call()
	dots <- as.list(match.call(expand.dots = FALSE)$...)
	data.null <- identical(data, parent.frame())
	data <- eval(mcall$data)
	xt <- substitute(x)
	xlab <- deparse(xt)
	if (!"xlab" %in% names(dots)) dots <- c(dots, list(xlab=xlab))
	x <- eval(mcall$x, data)
	yt <- substitute(y)
	ylab <- deparse(yt)
	if (!"ylab" %in% names(dots)) dots <- c(dots, list(ylab=ylab))
	y <- eval(mcall$y, data)
	xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
	id <- eval(mcall$id, data)
	if (!"pch" %in% names(dots)) dots <- c(dots, list(pch=46))
	if (is.null(id)) {
		idlab <- NULL
		do.call('plot', c(list(x, y), dots))
	}
	else {
		idlab <- deparse(mcall$id)
		id <- factor(eval(mcall$id, data))
		do.call('plot', c(list(x, y, type='n'), dots))
		if (!"col" %in% names(dots)) dots <- c(dots, list(col='gray'))
		if (!"type" %in% names(dots)) dots <- c(dots, list(type='o'))
		tt <- by(data, id, function (z) do.call('lines', c(list(eval(yt, z) ~ eval(xt, z)), dots)))
	}
	sel <- rep(FALSE, length(x))
	res <- integer(0)
	title('click on outliers in plot - then right-click to escape')
	while (sum(sel) < n) {
		ans <- identify(x[!sel], y[!sel], n=1, plot=FALSE)
		if (!length(ans)) break
		ans <- which(!sel)[ans]
		do.call('points', c(list(x[ans], y[ans]), par.out))
		if (!is.null(id)) {
			text(x[ans], y[ans], label=id[ans], adj=c(0.5, 0.1))
			idt <- id==id[ans]
			ox <- order(x[idt])
			do.call('lines', c(list(x[idt][ox], y[idt][ox]), par.out))
		}	
		sel[ans] <- TRUE
		res <- c(res, ans)
	}
	res <- res[order(res)]
	if (data.null) res 
	else list(rows=res, data=data[res, c(idlab, xlab, ylab)])
}
