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

context("edgeworth code")#FOLDUP
test_that("edgeworth trivial normal",{#FOLDUP
	xvals <- seq(-2,2,length.out=501)
	raw.moments <- c(0,1,0,3,0)
	raw.cumuls <- moment2cumulant(raw.moments)
	d1 <- dapx_edgeworth(xvals, raw.cumuls)
	d2 <- dnorm(xvals)
	# they should match:
	expect_true(max(abs(d1-d2)) < 1e-5)

	p1 <- papx_edgeworth(xvals, raw.cumuls)
	p2 <- pnorm(xvals)
	# they should match:
	expect_true(max(abs(p1-p2)) < 1e-5)

	p1 <- papx_edgeworth(xvals, raw.cumuls, lower.tail=FALSE)
	p2 <- pnorm(xvals, lower.tail=FALSE)
	# they should match:
	expect_true(max(abs(p1-p2)) < 1e-5)

	xvals <- seq(-6,6,length.out=501)
	mu <- 2
	sigma <- 3
	raw.moments <- c(2,13,62,475,3182)
	raw.cumuls <- moment2cumulant(raw.moments)
	d1 <- dapx_edgeworth(xvals, raw.cumuls)
	d2 <- dnorm(xvals,mean=mu,sd=sigma)
	# they should match:
	expect_true(max(abs(d1-d2)) < 1e-5)

	p1 <- papx_edgeworth(xvals, raw.cumuls)
	p2 <- pnorm(xvals,mean=mu,sd=sigma)
	# they should match:
	expect_true(max(abs(p1-p2)) < 1e-5)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("edgeworth trivial chisquare",{#FOLDUP
	chidf <- 30
	ords <- seq(1,9)
	raw.moments <- exp(ords * log(2) + lgamma((chidf/2) + ords) - lgamma(chidf/2))
	raw.cumuls <- moment2cumulant(raw.moments)
	xvals <- seq(0.3,30,length.out=801)

	# for a one-sided distribution, like the chi-square
	d1 <- dapx_edgeworth(xvals, raw.cumuls, support=c(0,Inf))
	d2 <- dchisq(xvals,df=chidf)
	# they should match:
	expect_true(max(abs(d1-d2)) < 1e-4)

	p1 <- papx_edgeworth(xvals, raw.cumuls, support=c(0,Inf))
	p2 <- pchisq(xvals,df=chidf)
	# they should match:
	expect_true(max(abs(p1-p2)) < 1e-4)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("edgeworth trivial beta",{#FOLDUP
	shape1 <- 20
	shape2 <- 40
	ords <- seq(0,9)
	raw.moments <- cumprod((shape1 + ords) / (shape1 + shape2 + ords))
	raw.cumuls <- moment2cumulant(raw.moments)
	xvals <- seq(0.001,0.999,length.out=801)

	# beta like
	d1 <- dapx_edgeworth(xvals, raw.cumuls, support=c(0,1))
	d2 <- dbeta(xvals,shape1=shape1,shape2=shape2)
	#plot(xvals,d1)
	#lines(xvals,d2,col='red')
	# they _should_ match:
	expect_true(max(abs(d1-d2)) < 1e-4)

	p1 <- papx_edgeworth(xvals, raw.cumuls, support=c(0,1))
	p2 <- pbeta(xvals,shape1=shape1,shape2=shape2)
	#plot(xvals,p1)
	#lines(xvals,p2,col='red')
	# they should match:
	expect_true(max(abs(p1-p2)) < 1e-5)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD

context("edgeworth NA checks")#FOLDUP
test_that("edgeworth check NA",{#FOLDUP
	chidf <- 50
	ords <- seq(1,9)
	raw.moments <- exp(ords * log(2) + lgamma((chidf/2) + ords) - lgamma(chidf/2))

	xvals <- seq(0.8/chidf,chidf*1.25,length.out=201)
	suppo <- c(0,Inf)
	#for (bases in c('normal','gamma')) {
		#for (lp in c(FALSE,TRUE)) {
			#d1 <- dapx_edgeworth(xvals, raw.moments, support=suppo, basis=bases, log=lp)
			#expect_true(!any(is.na(d1)))
			#for (lt in c(FALSE,TRUE)) {
				#p1 <- papx_edgeworth(xvals, raw.moments, support=suppo, basis=bases, log.p=lp, lower.tail=lt)
				#expect_true(!any(is.na(p1)))
			#}
		#}
	#}
	# someday edgeworth should support different bases
		for (lp in c(FALSE,TRUE)) {
			d1 <- dapx_edgeworth(xvals, raw.moments, support=suppo, log=lp)
			expect_true(!any(is.na(d1)))
			for (lt in c(FALSE,TRUE)) {
				p1 <- papx_edgeworth(xvals, raw.moments, support=suppo, log.p=lp, lower.tail=lt)
				expect_true(!any(is.na(p1)))
			}
		}

	shape1 <- 10
	shape2 <- 30
	ords <- seq(0,8)
	raw.moments <- cumprod((shape1 + ords) / (shape1 + shape2 + ords))
	xvals <- seq(0.001,0.999,length.out=801)
	suppo <- c(0,1)
	#for (bases in c('normal','gamma','beta')) {
		#for (lp in c(FALSE,TRUE)) {
			#d1 <- dapx_edgeworth(xvals, raw.moments, support=suppo, basis=bases, log=lp)
			#expect_true(!any(is.na(d1)))
			#for (lt in c(FALSE,TRUE)) {
				#p1 <- papx_edgeworth(xvals, raw.moments, support=suppo, basis=bases, log.p=lp, lower.tail=lt)
				#expect_true(!any(is.na(p1)))
			#}
		#}
	#}
	# someday edgeworth should support different bases
		for (lp in c(FALSE,TRUE)) {
			d1 <- dapx_edgeworth(xvals, raw.moments, support=suppo, log=lp)
			expect_true(!any(is.na(d1)))
			for (lt in c(FALSE,TRUE)) {
				p1 <- papx_edgeworth(xvals, raw.moments, support=suppo, log.p=lp, lower.tail=lt)
				expect_true(!any(is.na(p1)))
			}
		}

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD

context("edgeworth tail checks")#FOLDUP
test_that("edgeworth check tail",{#FOLDUP
	chidf <- 50
	ords <- seq(1,9)
	raw.moments <- exp(ords * log(2) + lgamma((chidf/2) + ords) - lgamma(chidf/2))

	xvals <- seq(0.8/chidf,chidf*1.25,length.out=201)
	suppo <- c(0,Inf)
	#for (bases in c('normal','gamma')) {
		#p1 <- papx_edgeworth(xvals, raw.moments, support=suppo, basis=bases, log.p=FALSE, lower.tail=TRUE)
		#p2 <- papx_edgeworth(xvals, raw.moments, support=suppo, basis=bases, log.p=FALSE, lower.tail=FALSE)
		#expect_true(max(abs(p1 + p2 - 1.0)) < 1e-5)
	#}
	# someday edgeworth should support different bases
	p1 <- papx_edgeworth(xvals, raw.moments, support=suppo, log.p=FALSE, lower.tail=TRUE)
	p2 <- papx_edgeworth(xvals, raw.moments, support=suppo, log.p=FALSE, lower.tail=FALSE)
	expect_true(max(abs(p1 + p2 - 1.0)) < 1e-5)

	shape1 <- 10
	shape2 <- 30
	ords <- seq(0,8)
	raw.moments <- cumprod((shape1 + ords) / (shape1 + shape2 + ords))
	xvals <- seq(0.001,0.999,length.out=801)
	suppo <- c(0,1)
	#for (bases in c('normal','gamma','beta')) {
		#p1 <- papx_edgeworth(xvals, raw.moments, support=suppo, basis=bases, log.p=FALSE, lower.tail=TRUE)
		#p2 <- papx_edgeworth(xvals, raw.moments, support=suppo, basis=bases, log.p=FALSE, lower.tail=FALSE)
		#expect_true(max(abs(p1 + p2 - 1.0)) < 1e-5)
	#}
	# someday edgeworth should support different bases
	p1 <- papx_edgeworth(xvals, raw.moments, support=suppo, log.p=FALSE, lower.tail=TRUE)
	p2 <- papx_edgeworth(xvals, raw.moments, support=suppo, log.p=FALSE, lower.tail=FALSE)
	expect_true(max(abs(p1 + p2 - 1.0)) < 1e-5)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
