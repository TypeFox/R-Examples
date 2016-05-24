`createWindow` <-
function(n) {
	ncol = trunc(sqrt(n))
	nrows = ceiling(n / ncol)
	par(mfrow = c(nrows, ncol))
}

