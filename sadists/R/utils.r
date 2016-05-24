# Copyright 2014-2014 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of sadists.
#
# sadists is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# sadists is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with sadists.  If not, see <http://www.gnu.org/licenses/>.

# Created: 2014.02.23
# Copyright: Steven E. Pav, 2014
# Author: Steven E. Pav
# Comments: Steven E. Pav

# compute the log (gamma(a) / gamma(b))
lgamrat <- function(a,b) {
	return(lgamma(a) - lgamma(b))
}
# compute the ratio gamma(a) / gamma(b)
gamrat <- function(a,b) {
	return(exp(lgamrat(a,b)))
}

# rchisq is currently broken for df=0,ncp=0
# this is a not-recycled chisquare generator
# which irons over the df=0,ncp=0 case.
unbroken_rchisq <- function(n,df,ncp=0) {
	retv <- rchisq(n,df=df,ncp=ncp)
	retv[df==0 & ncp==0] <- 0
	retv
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
