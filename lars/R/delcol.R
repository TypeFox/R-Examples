delcol <-
function(r, z, k = p)
{
	p <- dim(r)[1]
	r <- r[,  - k, drop = FALSE]
	z <- as.matrix(z)
	dz <- dim(z)
	storage.mode(r) <- storage.mode(z) <- "double"
	.Fortran("delcol",
		r,
		as.integer(p),
		as.integer(k),
		z,
		as.integer(dz[1]),
		as.integer(dz[2]),
                PACKAGE="lars")[c(1, 4)]
}

