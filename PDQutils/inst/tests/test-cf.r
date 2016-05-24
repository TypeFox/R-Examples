# Copyright 2015 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# This file is part of PDQutils.
#
# PDQutils is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PDQutils is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with PDQutils.  If not, see <http://www.gnu.org/licenses/>.

# env var:
# nb: 
# see also:
# todo:
# changelog: 
#
# Created: 2015.06.06
# Copyright: Steven E. Pav, 2015-2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}
THOROUGHNESS <- getOption('test.thoroughness',1.0)
#UNFOLD

context("CF code")#FOLDUP
test_that("CF trivial normal",{#FOLDUP
	pvals <- seq(0.001,0.999,length.out=501)
	q1 <- qapx_cf(pvals, c(0,1,0,0,0,0,0))
	q2 <- qnorm(pvals)

	# they should match:
	expect_true(max(abs(q1-q2)) < 1e-5)

	mu <- 2
	sigma <- 3
	raw.moments <- c(2,13,62,475,3182)
	raw.cumuls <- moment2cumulant(raw.moments)
	q1 <- qapx_cf(pvals, raw.cumuls)
	q2 <- qnorm(pvals,mean=mu,sd=sigma)
	# they should match:
	expect_true(max(abs(q1-q2)) < 1e-5)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("CF approximate chisquare",{#FOLDUP
	chidf <- 30
	ords <- seq(1,9)
	raw.moments <- exp(ords * log(2) + lgamma((chidf/2) + ords) - lgamma(chidf/2))
	raw.cumuls <- moment2cumulant(raw.moments)
	pvals <- seq(0.001,0.999,length.out=501)

	# for a one-sided distribution, like the chi-square
	q1 <- qapx_cf(pvals, raw.cumuls, support=c(0,Inf))
	q2 <- qchisq(pvals,df=chidf)
	# they should match:
	expect_true(max(abs(q1-q2)) < 1e-4)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("GC approximate beta",{#FOLDUP
	shape1 <- 20
	shape2 <- 40
	ords <- seq(0,6)
	raw.moments <- cumprod((shape1 + ords) / (shape1 + shape2 + ords))
	raw.cumuls <- moment2cumulant(raw.moments)
	pvals <- seq(0.001,0.999,length.out=501)

	# beta like
	q1 <- qapx_cf(pvals, raw.cumuls, support=c(0,1))
	q2 <- qbeta(pvals,shape1=shape1,shape2=shape2)
	# they _should_ match:
	expect_true(max(abs(q1-q2)) < 1e-4)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("GC rapx",{#FOLDUP
	set.seed(123120)
	r1 <- rapx_cf(50000, c(0,1,0,0,0,0,0))
	expect_true(abs(mean(r1)) < 1e-2)
	expect_true(abs(sd(r1) - 1) < 1e-2)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
