# /usr/bin/r
#
# Copyright 2015-2015 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of madness.
#
# madness is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# madness is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with madness.  If not, see <http://www.gnu.org/licenses/>.
#
# Created: 2015.11.18
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# internal utilities
.check_common_xtag <- function(e1,e2) {
	stopifnot(is.null(e1@xtag) || is.null(e2@xtag) || (e1@xtag == e2@xtag))
	if (is.null(e1@xtag)) {
		return(e2@xtag) # nocov
	} else {
		return(e1@xtag) # nocov
	}
}
.get_a_varx <- function(e1,e2) {
	if (is.null(e1@varx)) {
		return(e2@varx) # nocov
	} else {
		return(e1@varx) # nocov
	}
}

# apply commutator to Z, whose rows correspond to 
# vec(val)
# note that we require that tval has already been twiddled
# ('t') but Z has not.
.do_commutator <- function(tval,Z) {
	# we could use the commutation matrix from matrixcalc, 
	# but probably way faster to just reorder:
	newidx <- (row(tval) - 1) * ncol(tval) + col(tval)
	Z[newidx,,drop=FALSE] 
}

# quadratic equation roots, but beware roundoff.
.quadeq <- function(a,b,c) { # nocov start
	if (b < 0) {
		r1 <- (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
		r2 <- c / (a * r1)
	} else if (b > 0) {
		r1 <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
		r2 <- c / (a * r1)
	} else {
		r1 <- sqrt(c/a)
		r2 <- -r1
	}
	retv <- sort(c(r1,r2))
} # nocov end


#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
