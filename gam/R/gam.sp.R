"gam.sp" <-
function(x, knots, nknots, coef, smallest, scale)
{
	nas <- is.na(x)
	xs <- as.double((x[!nas] - smallest)/scale)
	bad.left <- xs < 0
	bad.right <- xs > 1
	good <- !(bad.left | bad.right)
	y <- xs
	if(any(good)) {
		junk <- .Fortran("bvalus",
			as.integer(sum(good)),
			knots,
			coef,
			as.integer(nknots),
			xs[good],
			s = double(sum(good)),
			as.integer(0),
                                 PACKAGE="gam")
		y[good] <- junk$s
	}
	if(any(!good)) {
		end.fit <- .Fortran("bvalus",
			as.integer(2),
			knots,
			coef,
			as.integer(nknots),
			as.double(c(0, 1)),
			s = double(2),
			as.integer(0),
                        PACKAGE="gam")$s
		end.slopes <- .Fortran("bvalus",
			as.integer(2),
			knots,
			coef,
			as.integer(nknots),
			as.double(c(0, 1)),
			s = double(2),
			as.integer(1),
                        PACKAGE="gam")$s
		if(any(bad.left))
			y[bad.left] <- end.fit[1] + end.slopes[1] * (xs[
				bad.left])
		if(any(bad.right))
			y[bad.right] <- end.fit[2] + end.slopes[2] * (xs[
				bad.right] - 1)
	}
	pred <- x * 0
	pred[!nas] <- y
	pred
}
