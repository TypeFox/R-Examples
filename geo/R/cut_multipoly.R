#' Intersect or take complement (?) of polygons
#' 
#' Intersect or take complement (?) of one or many polygons with a single one.
#' 
#' 
#' @param x Polygon or polygons, seperated with NAs, hence 'multipoly'
#' @param xb Polygon to intersect with/complement from
#' @param in.or.out Whether to take the intersect (0) or complement of x in xb
#' (1). Default 0
#' @return List with compents: \item{x, y}{with coordinates of intersection or
#' complement (?)}
#' @note Check use of \code{in.or.out=1}, when there are many polygons in
#' \code{x}, how is the complement with \code{xb} taken? Needs elaboration.
#' Possibly drop in.or.out from argument list and fix it at 0 in call to
#' findcut.
#' @seealso Called by \code{\link{geopolygon}}, calls \code{\link{findcut}} and
#' \code{\link{prepare.line}}
#' @keywords manip
#' @export cut_multipoly
cut_multipoly <-
function(x, xb, in.or.out = 0)
{
	ind <- x$x[is.na(x$x)]
	if(length(ind) == 0) {
		x2 <- findcut(x, xb, in.or.out)
	}
	else {
		x2 <- list(x = NA, y = NA)
		ind <- prepare.line(x$x)
		for(i in 1:ind$nlx) {
			x1 <- list(x = x$x[ind$lx1[i]:ind$lx2[i]], y = x$y[
				ind$lx1[i]:ind$lx2[i]])
			x1 <- findcut(x1, xb, in.or.out)
			x2$x <- c(x2$x, NA, x1$x)
			x2$y <- c(x2$y, NA, x1$y)
		}
		x2$x <- x2$x[ - c(1:2)]
		x2$y <- x2$y[ - c(1:2)]
	}
	return(x2)
}

