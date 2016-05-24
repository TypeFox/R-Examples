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

# Created: 2015.02.08
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# we need a vectorized confluent hypergeometric function#FOLDUP
F11 <- Vectorize(hypergeo::genhypergeo,vectorize.args=c("U","L","z"))
ReF11 <- function(...) {
	return(Re(F11(...)))
}
#UNFOLD

# some moments of standard distributions:#FOLDUP

# compute the 1 through order.max raw, uncentered moment
# of the normal distribution with given mean and standard
# deviation
norm_moms <- function(mu=0,sigma=1,order.max=3) {
	retval <- rep(1,order.max)
	hermi <- orthopolynom::hermite.he.polynomials(order.max, normalized=FALSE)
	for (iii in c(1:order.max)) {
		cvals <- abs(coefficients(hermi[[iii+1]]))
		lvals <- mu^(0:iii) * sigma^(iii - (0:iii))
		retval[iii] <- sum(cvals * lvals)
	}
	return(retval)
}

# compute the 1 through order.max raw, uncentered moment
# of the (non-central) chisquare to power distribution.
chipow_moms <- function(df,ncp=0,pow=1,order.max=3,orders=1:order.max,log=FALSE) {
	retval <- chisq_moms(df=df,ncp=ncp,orders=pow * orders,log=log)
	return(retval)
}

# compute the 1 through order.max raw, uncentered moment
# of the (non-central) chi-square distribution with df d.f.
# c.f. Paolella 10.9
chisq_moms <- function(df,ncp=0,order.max=3,orders=1:order.max,log=FALSE) {
	hadf <- df/2.0
	hancp <- ncp/2.0
	if (log) {
		retval <- orders * log(2.0) + lgamrat(orders + hadf,hadf)
		if (ncp != 0) {
			retval <- retval - hancp + log(ReF11(orders + hadf,hadf,hancp))
		}
	} else {
		retval <- (2.0^(orders)) * gamrat(orders + hadf,hadf)
		if (ncp != 0) {
			retval <- (retval/exp(hancp)) * ReF11(orders + hadf,hadf,hancp)
		}
	}
	return(retval)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
