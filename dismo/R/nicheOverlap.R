# based on function Istat in package SDMTools by Jeremy VanDerWal
# adapted by Robert J. Hijmans 
# Date : October 2012
# Version 1.0
# Licence GPL v3


nicheOverlap <- function (x, y, stat='I', mask=TRUE, checkNegatives=TRUE) {
	s <- stack(x, y)

	stopifnot(stat %in% c('I', 'D'))
	# to assure that both have the same NA cells
	if (mask) {
		s <- mask(s,  x * y)
	}
	
	if (checkNegatives) {
		minv <- cellStats(s, 'min', na.rm=TRUE)
		if (any(minv < 0)) {
			stop('values of "x" and "y" should be non-negative')
		}
	}
	
	cs <- cellStats(s, 'sum', na.rm=TRUE)
	if (stat == 'I') {
		r <- overlay(s, fun=function(i,j) (sqrt(i / cs[1]) - sqrt(j / cs[2] ))^2)
		H2 <- cellStats(r, 'sum', na.rm=TRUE)
		1 - 0.5 * H2
	} else {
		r <- overlay(s, fun=function(i,j) abs((i / cs[1]) - (j / cs[2])))
		D <- cellStats(r, 'sum', na.rm=TRUE)
		1 - 0.5 * D
	
	}
}


