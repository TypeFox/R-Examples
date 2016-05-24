# Copyright 2014 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# This file is part of MarkowitzR.
#
# MarkowitzR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MarkowitzR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MarkowitzR.  If not, see <http://www.gnu.org/licenses/>.

# env var:
# nb: 
# see also:
# todo:
# changelog: 
#
# Created: 2014.05.19
# Copyright: Steven E. Pav, 2014-2014
# Author: Steven E. Pav
# Comments: Steven E. Pav

# just test if everything runs...

# helpers#FOLDUP
set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}
THOROUGHNESS <- getOption('test.thoroughness',1.0)
#UNFOLD

context("test API")#FOLDUP
test_that("mp_vcov runs",{#FOLDUP
	ngen <- ceiling(THOROUGHNESS * 32)
	alpha.floor = 0.001 + 0.003 * (THOROUGHNESS / (1 + THOROUGHNESS))

	set.char.seed("0b144107-4de8-4e00-95f7-d746db3aef8e")
	ope <- 253
	for (nyr in c(2,4)) {
		nday <- ceiling(ope * nyr)
		for (nstock in c(4,6)) {
			X <- matrix(rnorm(nday * nstock),ncol=nstock)
			expect_that(asym <- MarkowitzR::mp_vcov(X,fit.intercept=TRUE),
									not(throws_error()))
			for (psiz in c(1,2)) {
				Amat <- matrix(rnorm(1 * nstock),ncol=nstock)
				expect_that(asym <- MarkowitzR::mp_vcov(X,Jmat=Amat,
																								fit.intercept=TRUE),
										not(throws_error()))
				expect_that(asym <- MarkowitzR::mp_vcov(X,Gmat=Amat,
																								fit.intercept=TRUE),
										not(throws_error()))

			}
		}
	}
})#UNFOLD

#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:

