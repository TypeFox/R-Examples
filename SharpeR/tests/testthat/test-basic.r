# Copyright 2016 Steven E. Pav. All Rights Reserved.
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
# Created: 2016.03.12
# Copyright: Steven E. Pav, 2013-2016
# Author: Steven E. Pav
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}
THOROUGHNESS <- getOption('test.thoroughness',1.0)
#UNFOLD

context("code runs at all")#FOLDUP
test_that("sr_basic",{#FOLDUP
	set.char.seed("653c0960-a65c-4091-b916-cb9bfc21a21b")
	X <- matrix(rnorm(1000*1),ncol=1)
	Y <- matrix(rnorm(length(X)),ncol=1)
	asrX <- as.sr(X)
	asrY <- as.sr(Y)
	print(summary(asrX))
	print(summary(asrY))

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("sropt_test",{#FOLDUP
	set.char.seed("1f9001aa-7457-4492-b4a0-ce7b8ae55049")
	X <- matrix(rnorm(1000*10),ncol=10)
	sroX <- as.sropt(X)
	print(summary(sroX))

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
