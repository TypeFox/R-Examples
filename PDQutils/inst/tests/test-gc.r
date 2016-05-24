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

context("GC code runs")#FOLDUP
test_that("GC trivial normal",{#FOLDUP
	xvals <- seq(-2,2,length.out=501)
	raw.moments <- c(0,1,0,3,0)
	d1 <- dapx_gca(xvals, raw.moments, basis='normal')
	d2 <- dnorm(xvals)
	# they should match:
	expect_true(max(abs(d1-d2)) < 1e-5)

	p1 <- papx_gca(xvals, raw.moments)
	p2 <- pnorm(xvals)
	# they should match:
	expect_true(max(abs(p1-p2)) < 1e-5)

	xvals <- seq(-6,6,length.out=501)
	mu <- 2
	sigma <- 3
	raw.moments <- c(2,13,62,475,3182)
	d1 <- dapx_gca(xvals, raw.moments, basis='normal')
	d2 <- dnorm(xvals,mean=mu,sd=sigma)
	# they should match:
	expect_true(max(abs(d1-d2)) < 1e-5)

	p1 <- papx_gca(xvals, raw.moments, basis='normal')
	p2 <- pnorm(xvals,mean=mu,sd=sigma)
	# they should match:
	expect_true(max(abs(p1-p2)) < 1e-5)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("GC trivial chisquare",{#FOLDUP
	chidf <- 30
	ords <- seq(1,9)
	raw.moments <- exp(ords * log(2) + lgamma((chidf/2) + ords) - lgamma(chidf/2))
	xvals <- seq(0.3,30,length.out=801)

	# for a one-sided distribution, like the chi-square
	d1 <- dapx_gca(xvals, raw.moments, support=c(0,Inf), basis='gamma')
	d2 <- dchisq(xvals,df=chidf)
	# they should match:
	expect_true(max(abs(d1-d2)) < 1e-5)

	p1 <- papx_gca(xvals, raw.moments, support=c(0,Inf), basis='gamma')
	p2 <- pchisq(xvals,df=chidf)
	# they should match:
	expect_true(max(abs(p1-p2)) < 1e-5)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("GC trivial beta",{#FOLDUP
	shape1 <- 20
	shape2 <- 40
	ords <- seq(0,6)
	raw.moments <- cumprod((shape1 + ords) / (shape1 + shape2 + ords))
	xvals <- seq(0.001,0.999,length.out=801)

	# beta like
	#d1 <- dapx_gca(xvals, raw.moments, support=c(0,1), basis='beta', basepar=list(shape1=shape1,shape2=shape2))
	d1 <- dapx_gca(xvals, raw.moments, support=c(0,1), basis='beta')
	d2 <- dbeta(xvals,shape1=shape1,shape2=shape2)
	#plot(xvals,d1)
	#lines(xvals,d2,col='red')
	# they _should_ match:
	expect_true(max(abs(d1-d2)) < 1e-5)

	p1 <- papx_gca(xvals, raw.moments, support=c(0,1), basis='beta')
	p2 <- pbeta(xvals,shape1=shape1,shape2=shape2)
	#plot(xvals,p1)
	#lines(xvals,p2,col='red')
	# they should match:
	expect_true(max(abs(p1-p2)) < 1e-5)

	# sentinel
	expect_true(TRUE)
})#UNFOLD

test_that("GC check runs",{#FOLDUP

	# these are nfg for now....
	xvals <- seq(-1,1,length.out=501)
	raw.moments <- c(0,1,0,2,0)
	d1 <- dapx_gca(xvals, raw.moments, basis='wigner')
	p1 <- papx_gca(xvals, raw.moments, basis='wigner')

	xvals <- seq(-1,1,length.out=501)
	raw.moments <- c(0,1,3,4,0)
	d1 <- dapx_gca(xvals, raw.moments, basis='arcsine')
	p1 <- papx_gca(xvals, raw.moments, basis='arcsine')

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD

context("GC NA checks")#FOLDUP
test_that("GC check NA",{#FOLDUP
	chidf <- 50
	ords <- seq(1,9)
	raw.moments <- exp(ords * log(2) + lgamma((chidf/2) + ords) - lgamma(chidf/2))

	xvals <- seq(0.8/chidf,chidf*1.25,length.out=201)
	suppo <- c(0,Inf)
	for (bases in c('normal','gamma')) {
		for (lp in c(FALSE,TRUE)) {
			d1 <- dapx_gca(xvals, raw.moments, support=suppo, basis=bases, log=lp)
			expect_true(!any(is.na(d1)))
			for (lt in c(FALSE,TRUE)) {
				p1 <- papx_gca(xvals, raw.moments, support=suppo, basis=bases, log.p=lp, lower.tail=lt)
				expect_true(!any(is.na(p1)))
			}
		}
	}

	shape1 <- 10
	shape2 <- 30
	ords <- seq(0,8)
	raw.moments <- cumprod((shape1 + ords) / (shape1 + shape2 + ords))
	xvals <- seq(0.001,0.999,length.out=801)
	suppo <- c(0,1)
	for (bases in c('normal','gamma','beta')) {
		for (lp in c(FALSE,TRUE)) {
			d1 <- dapx_gca(xvals, raw.moments, support=suppo, basis=bases, log=lp)
			expect_true(!any(is.na(d1)))
			for (lt in c(FALSE,TRUE)) {
				p1 <- papx_gca(xvals, raw.moments, support=suppo, basis=bases, log.p=lp, lower.tail=lt)
				expect_true(!any(is.na(p1)))
			}
		}
	}

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD

context("GC tail checks")#FOLDUP
test_that("GC check tail",{#FOLDUP
	chidf <- 50
	ords <- seq(1,9)
	raw.moments <- exp(ords * log(2) + lgamma((chidf/2) + ords) - lgamma(chidf/2))

	xvals <- seq(0.8/chidf,chidf*1.25,length.out=201)
	suppo <- c(0,Inf)
	for (bases in c('normal','gamma')) {
		p1 <- papx_gca(xvals, raw.moments, support=suppo, basis=bases, log.p=FALSE, lower.tail=TRUE)
		p2 <- papx_gca(xvals, raw.moments, support=suppo, basis=bases, log.p=FALSE, lower.tail=FALSE)
		expect_true(max(abs(p1 + p2 - 1.0)) < 1e-5)
	}

	shape1 <- 10
	shape2 <- 30
	ords <- seq(0,8)
	raw.moments <- cumprod((shape1 + ords) / (shape1 + shape2 + ords))
	xvals <- seq(0.001,0.999,length.out=801)
	suppo <- c(0,1)
	for (bases in c('normal','gamma','beta')) {
		p1 <- papx_gca(xvals, raw.moments, support=suppo, basis=bases, log.p=FALSE, lower.tail=TRUE)
		p2 <- papx_gca(xvals, raw.moments, support=suppo, basis=bases, log.p=FALSE, lower.tail=FALSE)
		expect_true(max(abs(p1 + p2 - 1.0)) < 1e-5)
	}

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD

context("moments code sane")#FOLDUP
test_that("moments basics",{#FOLDUP
	raw.mom1 <- c(0,1,0,3,0)
	raw.cumul1 <- moment2cumulant(raw.mom1)
	raw.mom2 <- cumulant2moment(raw.cumul1)
	raw.cumul2 <- moment2cumulant(raw.mom2)

	expect_true(max(abs(raw.mom1 - raw.mom2)) < 1e-9)
	expect_true(max(abs(raw.cumul1 - raw.cumul2)) < 1e-9)

	chidf <- 30
	ords <- seq(1,9)
	raw.mom1 <- exp(ords * log(2) + lgamma((chidf/2) + ords) - lgamma(chidf/2))
	raw.cumul1 <- moment2cumulant(raw.mom1)
	raw.mom2 <- cumulant2moment(raw.cumul1)
	raw.cumul2 <- moment2cumulant(raw.mom2)

	expect_true(max(abs(raw.mom1 - raw.mom2)) < 1e-9)
	expect_true(max(abs(raw.cumul1 - raw.cumul2)) < 1e-9)

	shape1 <- 20
	shape2 <- 40
	ords <- seq(0,6)
	raw.mom1 <- cumprod((shape1 + ords) / (shape1 + shape2 + ords))
	raw.cumul1 <- moment2cumulant(raw.mom1)
	raw.mom2 <- cumulant2moment(raw.cumul1)
	raw.cumul2 <- moment2cumulant(raw.mom2)

	expect_true(max(abs(raw.mom1 - raw.mom2)) < 1e-9)
	expect_true(max(abs(raw.cumul1 - raw.cumul2)) < 1e-9)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
