##
##  v e c t o r f i e l d . R
##


vectorfield <- function(fun, xlim, ylim, n = 16,
                        scale = 0.05, col = "green", ...) {
    stopifnot(is.numeric(xlim), length(xlim) == 2,
              is.numeric(ylim), length(ylim) == 2)

	xpts <- linspace(xlim[1], xlim[2], n)
	ypts <- linspace(ylim[1], ylim[2], n)

	M <- meshgrid(xpts, ypts)
	x <- M$X
	y <- M$Y

	px = matrix(1, nrow=n , ncol=n)
	py = fun(x, y);

	plot(xlim, ylim, type="n"); grid()
	quiver(x, y, px, py, scale = scale, col = col, ...)
}


quiver <- function(x, y, u, v,
                    scale = 0.05, angle = 10, length = 0.1, ...) {
    stopifnot(is.numeric(x), is.numeric(y), is.numeric(u), is.numeric(v))
    

	arrows(x, y, x+scale*u, y+scale*v, angle=10, length=length, ...)
}
