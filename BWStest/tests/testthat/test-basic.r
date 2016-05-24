# Copyright 2016 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# This file is part of BWStest.
#
# BWStest is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BWStest is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BWStest.  If not, see <http://www.gnu.org/licenses/>.

# env var:
# nb: 
# see also:
# todo:
# changelog: 
#
# Created: 2016.04.06
# Copyright: Steven E. Pav, 2016-2016
# Author: Steven E. Pav
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}
THOROUGHNESS <- getOption('test.thoroughness',1.0)
#UNFOLD

context("code runs at all")#FOLDUP
test_that("runs without error",{#FOLDUP
	set.char.seed("46316ac4-f424-43f9-b5c2-1f1897e5abbd")
	x <- rnorm(100)
	y <- rnorm(100)

	b <- bws_stat(x,y)
	b <- bws_stat(x,rnorm(1))
	b <- bws_stat(rnorm(1),y)

	bvals <- replicate(500,bws_stat(rnorm(50),rnorm(40)))
	pvals <- bws_cdf(bvals)

	htest <- bws_test(x,y)

	expect_warning(bws_test(rnorm(5),rnorm(5)))
	expect_warning(bws_test(rnorm(50),rnorm(5)))

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("murakami computations",{#FOLDUP
	set.char.seed("f1a7627e-3545-4fdf-a5bf-10cf9a4bb396")
	for (flavor in 0:5) {
		statv <- murakami_stat(rnorm(10),rnorm(10),flavor)
		allev <- murakami_stat_perms(3,3,flavor)
		allev <- murakami_stat_perms(6,3,flavor)

		# hope this errors, or it will consume all your memory!
		expect_error(allev <- murakami_stat_perms(40,40,flavor))
	}
	expect_error(statv <- murakami_stat(rnorm(10),rnorm(10),flavor=-1))
	expect_error(statv <- murakami_stat(rnorm(10),rnorm(10),flavor=6))

	# all for coverage
	for (lower_tail in c(FALSE,TRUE)) {
		for (flavor in c(0:5)) {
			cdfv <- murakami_cdf(seq(-1,1,length.out=101),n1=5,n2=5,flavor=flavor,lower_tail=lower_tail)
			cdfv <- murakami_cdf(seq(-1,1,length.out=101),n1=5,n2=4,flavor=flavor,lower_tail=lower_tail)
			cdfv <- murakami_cdf(seq(-1,1,length.out=101),n1=4,n2=5,flavor=flavor,lower_tail=lower_tail)
		}
		expect_error(cdfv <- murakami_cdf(seq(-1,1,length.out=101),n1=-5,n2=5,flavor=0L,lower_tail=lower_tail))
		expect_error(cdfv <- murakami_cdf(seq(-1,1,length.out=101),n1=4,n2=5,flavor=-1,lower_tail=lower_tail))
		expect_error(cdfv <- murakami_cdf(seq(-1,1,length.out=101),n1=4,n2=5,flavor=6,lower_tail=lower_tail))
	}

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("bws_test",{#FOLDUP
	set.char.seed("140d1002-1cc4-426a-8489-d7cd62988d0b")

	x <- rnorm(8)
	y <- rnorm(8)

	for (mth in c('default','BWS','Neuhauser','B1','B2','B3','B4','B5')) {
		hval <- bws_test(x,y,method=mth,alternative='two.sided')
	}

	for (mth in c('Neuhauser','B2')) {
		for (alt in c("two.sided","greater","less")) {
			hval <- bws_test(x,y,method=mth,alternative=alt)
		}
	}
	for (mth in c('BWS','B1','B3','B4','B5')) {
		for (alt in c("greater","less")) {
			expect_error(hval <- bws_test(x,y,method=mth,alternative=alt))
		}
	}

	# better get warnings out of this due to small sample size...
	x <- rnorm(20)
	y <- rnorm(20)

	for (mth in c('B3','B4','B5')) {
		for (alt in c("two.sided")) {
			expect_warning(hval <- bws_test(x,y,method=mth,alternative=alt))
		}
	}

	# sentinel
	expect_true(TRUE)
})#UNFOLD

# 2FIX: check the effects of NA
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
