#' Find intersection or complement
#' 
#' Find intersection or compliment of two polygons.
#' 
#' 
#' @param x Polygon
#' @param xb Polygon to intersect with/complement from
#' @param in.or.out Whether to take intersect of \code{x} and \code{xb} (0) or
#' complement of \code{x} in \code{xb} (1). Default 0.
#' @return Returns a list of \item{x, y}{Coordinate of intersect or compliment}
#' \item{nxr}{Number/index of returned coordinates in \code{xb} (?)}
#' @note Needs elaboration.
#' @seealso Called by \code{\link{cut_multipoly}}, \code{\link{geointersect}}
#' and \code{\link{reitaplott}}; calls \code{\link{find.hnit}} and
#' \code{\link{geoinside}}.
#' @keywords manip logic
#' @export findcut
findcut <-
function(x, xb, in.or.out)
{
	if(!is.data.frame(x))
		x <- data.frame(x = x$x, y = x$y)
	xr <- yr <- mark <- side <- s <- t <- rep(0, (length(x$y) + length(
		xb$y)))
	nxr <- 0
	ab <- ab1 <- rep(0, length(xb$x))
	xr <- .C("define_poly", PACKAGE = "geo", 
		as.double(x$x),
		as.double(x$y),
		as.double(xb$x),
		as.double(xb$y),
		as.double(xr),
		as.double(yr),
		as.integer(length(x$y)),
		as.integer(length(xb$y)),
		as.integer(nxr),
		as.integer(mark),
		as.integer(side),
		as.double(s),
		as.double(t),
		as.double(ab),
		as.double(ab1),
		as.integer(in.or.out))
	nxr <- xr[[9]]
	yr <- xr[[6]][1:nxr]
	mark <- xr[[10]][1:nxr]
	side <- xr[[11]][1:nxr]
	s <- xr[[12]][1:nxr]
	t <- xr[[12]][1:nxr]
	xr <- xr[[5]][1:nxr]
	ind <- c(1:nxr)
	ind2 <- ind[mark == 2]
	i <- geoinside(x[1,  ], reg = xb, option = 3, col.names = c("x", "y"))
	if(in.or.out == 1)
		i <- !i
	if(i == 1 && length(ind2) == 0)
		return(list(x = x$x, y = x$y))
	if(i == 0 && length(ind2) == 0)
		return(invisible())
	if(ind2[1] == 1)
		ind1 <- c(1:nxr)
	else ind1 <- c(ind2[1]:nxr, 1:(ind2[1] - 1))
	xr <- xr[ind1]
	yr <- yr[ind1]
	mark <- mark[ind1]
	side <- side[ind1]
	s <- s[ind1]
	t <- t[ind1]
	h1 <- side + s + 1
	inn <- ifelse(mark == 2, 1, 0)
	ind1 <- ind[mark == 1 | mark == 2]
	nr <- ind[mark == 1 | mark == 2]
	h <- h1[ind1]
	n <- length(h)
	if(n < 2)
		return(invisible())
	# vidbot i profun
	s <- matrix(0, n, 3)
	s[, 2] <- match(sort(h), h)
	s[, 1] <- c(s[n, 2], s[1:(n - 1), 2])
	s[, 3] <- c(s[2:n, 2], s[1, 2])
	o <- match(h, sort(h))
	s <- s[o,  ]
	up <- rep(0, nrow(s))
	pt <- h[s[, 2]] - 0.0001
	i <- geoinside(find.hnit(pt, xb), reg = x, option = 0, col.names = c(
		"x", "y"))
	if(in.or.out == 1) {
		i1 <- c(1:nrow(find.hnit(pt, xb)))
		i <- i1[is.na(match(i1, i))]
	}
	if(length(i) > 0) {
		s[ - i, 1] <- s[ - i, 3]
		up[ - i] <- 1
	}
	s <- matrix(c(s[, 2], s[, 1]),  , 2)
	s1 <- matrix(0, length(h1), 2)
	s1[ind1,  ] <- s
	up1 <- buid <- rep(0, length(h1))
	up1[ind1] <- up
	s1[, 2] <- match(s1[, 2], s1[, 1])
	s1[, 1] <- 1:nrow(s1)
	nxr <- 0
	xr1 <- yr1 <- rep(0, (length(x$x) + length(xb$x)))
	x <- .C("post_filter", PACKAGE = "geo", 
		as.integer(s1[, 2]),
		as.integer(side),
		as.integer(up1),
		as.integer(mark),
		as.double(xr),
		as.double(yr),
		as.integer(buid),
		as.integer(nrow(s1)),
		as.double(xb$x),
		as.double(xb$y),
		as.integer(length(xb$y)),
		as.double(xr1),
		as.double(yr1),
		as.integer(nxr))
	nxr <- x[[14]]
	xr <- x[[12]][1:nxr]
	yr <- x[[13]][1:nxr]
	ind <- c(1:nxr)
	ind <- ind[xr < -999998]
	if(length(ind) > 0)
		xr[ind] <- yr[ind] <- NA
	return(list(x = xr, y = yr, nxr = nxr))
}

