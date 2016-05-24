# Copyright 2014-2015 Steven E. Pav. All Rights Reserved.
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

# Created: 2015.02.23
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# compute the 1 through order.max raw, cumulants 
# of the normal distribution with given mean and sd.
norm_cumuls <- function(mean=0,sd=1,order.max=3) {
	retval <- rep(0,order.max)
	if (order.max > 0)
		retval[1] <- mean
	if (order.max > 1)
		retval[2] <- sd
	return(retval)
}

# compute the 1 through order.max raw cumulants of
# the (non-central) chi-square to the pow power distribution.
chipow_cumuls <- function(df,ncp=0,pow=1,order.max=3) {
	retval <- moment2cumulant(chipow_moms(df,ncp=ncp,pow=pow,order.max=order.max))
	return(retval)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
