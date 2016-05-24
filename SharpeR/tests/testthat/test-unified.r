# Copyright 2013 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# This file is part of SharpeR.
#
# SharpeR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SharpeR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SharpeR.  If not, see <http://www.gnu.org/licenses/>.

# env var:
# nb: 
# see also:
# todo:
# changelog: 
#
# Created: 2013.12.23
# Copyright: Steven E. Pav, 2013-2013
# Author: Steven E. Pav
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}
THOROUGHNESS <- getOption('test.thoroughness',1.0)
#UNFOLD

context("test unified Markowitz")#FOLDUP
test_that("sm/ism_vcov runs",{#FOLDUP
	set.char.seed("1023fe77-6323-478f-bea7-304d8008dde8")

	X <- matrix(rnorm(1000*4),ncol=4)
	# just want to make sure these all run:
	foo <- sm_vcov(X)
	baz <- ism_vcov(X)
	foo <- sm_vcov(X,fit.intercept=FALSE)
	baz <- ism_vcov(X,fit.intercept=FALSE)
	foo <- sm_vcov(X,vcov.func="normal")
	baz <- ism_vcov(X,vcov.func="normal")

	expect_error(sm_vcov(X,vcov.func="fnordmal"))
	expect_error(sm_vcov(X,vcov.func="normal",fit.intercept=FALSE))

	# sentinel
	expect_true(TRUE)
})#UNFOLD

#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
