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
# Created: 2013.04.26
# Copyright: Steven E. Pav, 2013-2013
# Author: Steven E. Pav
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}
THOROUGHNESS <- getOption('test.thoroughness',1.0)
#UNFOLD

context("code runs at all")#FOLDUP
test_that("basic sr functionality",{
	set.char.seed("50e02c82-fdfc-4423-93ef-cbd350a65391")
	x <- matrix(rnorm(100*2),ncol=4)
	mysr <- as.sr(x)
	print(mysr)
	print(format(mysr))
	dummy <- reannualize(mysr,new.ope=2)
	expect_true(is.sr(mysr))
	myse <- se(mysr,"t")
	myse <- se(mysr,"Lo")

	ci1 <- confint(mysr)
	ci2 <- confint(mysr, level=0.9)
	ci3 <- confint(mysr, level.lo=0.1, level.hi=0.95)
	ci4 <- confint(mysr, type="exact")
	ci5 <- confint(mysr, type="t")
	ci6 <- confint(mysr, type="Z")

	# via data.frame
	xdf <- data.frame(a=rnorm(500,1),b=rnorm(500,1));
	mysr <- as.sr(xdf)
	print(mysr)
	print(format(mysr))
	expect_true(is.sr(mysr))

	# via lm
	xdf <- data.frame(AAPL=rnorm(500,1),SPY=rnorm(500,1));
	mylm <- lm(AAPL ~ SPY,data=xdf)
	mysr <- as.sr(mylm)
	print(mysr)
	print(format(mysr))
	expect_true(is.sr(mysr))

	# via xts
	if (require(xts)) {
		xxts <- xts(x,order.by=as.Date(1:(dim(x)[1]),origin="1990-01-01"))
		mysr <- as.sr(xxts)
		print(mysr)
		print(format(mysr))
		expect_true(is.sr(mysr))
		dummy <- reannualize(mysr,new.ope=2)

		if (require(timeSeries)) {
			ats <- as.timeSeries(xxts)
			mysr <- as.sr(ats)
			print(mysr)
			print(format(mysr))
			expect_true(is.sr(mysr))
			dummy <- reannualize(mysr,new.ope=2)
		}
	}

	expect_true(TRUE)
})

test_that("basic sropt functionality",{
	set.char.seed("943c439d-8f66-4f5e-a2a5-a9af7cf1b787")
	x <- matrix(rnorm(1000*10),ncol=10)
	mysr <- as.sropt(x)
	print(mysr)
	print(format(mysr))
	dummy <- reannualize(mysr,new.ope=2)
	expect_true(is.sropt(mysr))
	ci1 <- confint(mysr)
	ci2 <- confint(mysr, level=0.9)
	ci3 <- confint(mysr, level.lo=0.1, level.hi=0.95)

	z1 <- inference(mysr, type="KRS")
	z2 <- inference(mysr, type="MLE")
	z3 <- inference(mysr, type="unbiased")

	sricv <- sric(mysr)
	sricv <- sric(dummy)

	# via xts
	if (require(xts)) {
		xxts <- xts(x,order.by=as.Date(1:(dim(x)[1]),origin="1990-01-01"))
		mysr <- as.sropt(xxts)
		print(mysr)
		print(format(mysr))
		expect_true(is.sropt(mysr))
		dummy <- reannualize(mysr,new.ope=2)
	}

	expect_true(TRUE)
})
test_that("basic del_sropt functionality",{
	set.char.seed("cdbfbd76-29f4-454f-b8a1-e3502ba287d7")
	x <- matrix(rnorm(1000*10,mean=3e-4),ncol=10)
	G <- matrix(rnorm(2*10),ncol=10)
	mysr <- as.del_sropt(x,G,drag=0,ope=1,epoch="yr")
	print(mysr)
	print(format(mysr))
	# not yet:
	#dummy <- reannualize(mysr,new.ope=2)
	expect_true(is.del_sropt(mysr))
	ci1 <- confint(mysr)
	ci2 <- confint(mysr, level=0.9)
	ci3 <- confint(mysr, level.lo=0.1, level.hi=0.95)

	z1 <- inference(mysr, type="KRS")
	z2 <- inference(mysr, type="MLE")
	z3 <- inference(mysr, type="unbiased")

	# via xts
	if (require(xts)) {
		xxts <- xts(x,order.by=as.Date(1:(dim(x)[1]),origin="1990-01-01"))
		mysr <- as.del_sropt(xxts,G,drag=0,ope=1,epoch="yr")
		print(mysr)
		print(format(mysr))
		expect_true(is.del_sropt(mysr))
		#dummy <- reannualize(mysr,new.ope=2)
	}

	expect_true(TRUE)
})
test_that("basic sr_vcov functionality",{
	set.char.seed("de096679-85cc-438b-8335-96a9940a9021")
	for (p in c(1:3)) {
		X <- matrix(rnorm(1000*p,mean=3e-4),ncol=p)
		S <- sr_vcov(X)
	}
	X <- rnorm(1000)
	S <- sr_vcov(X)

	# sentinel:
	expect_true(TRUE)
})
#UNFOLD
context("estimation functions: confint coverage")#FOLDUP
test_that("confint.sr coverage",{#FOLDUP
	set.char.seed("066dfa96-6dd7-4d14-ab74-49e81b3afd83")

	ngen <- ceiling(THOROUGHNESS * 32)
	alpha.tol = 0.05 + 0.10 / THOROUGHNESS

	ope <- 253
	for (nyr in c(3,6,9)) {
		df <- ceiling(ope * nyr)
		for (type in c("exact","t","Z")) {
			for (zeta in c(-1.0,0.0,1.0,2.0)) {
				rvs <- rsr(ngen, df, zeta, ope)
				roll.own <- sr(sr=rvs,df=df,c0=0,ope=ope)
				for (nominal.coverage in c(0.90,0.95)) {
					aci <- confint(roll.own,level=nominal.coverage,type=type)
					coverage <- 1 - mean((zeta < aci[,1]) | (aci[,2] < zeta))
					expect_equal(coverage,nominal.coverage,tolerance=alpha.tol,scale=1)
				}
			}
		}
	}
})#UNFOLD
test_that("confint.sropt coverage",{#FOLDUP
	set.char.seed("50c7aa74-9cec-4fef-980e-c25bfede8260")

	ngen <- ceiling(THOROUGHNESS * 64)
	alpha.tol = 0.05 + 0.10 / THOROUGHNESS

	ope <- 253
	for (df1 in c(2,4,8,16)) {
		for (nyr in c(3,6,9)) {
			df2 <- ceiling(ope * nyr)
			for (zeta.s in c(0.0,1.0,3.0)) {
				rvs <- rsropt(ngen, df1, df2, zeta.s, ope)
				roll.own <- sropt(z.s=rvs,df1,df2,drag=0,ope=ope)
				for (nominal.coverage in c(0.90,0.95)) {
					aci <- confint(roll.own,level=nominal.coverage)
					coverage <- 1 - mean((zeta.s < aci[,1]) | (aci[,2] < zeta.s))
					#cat(sprintf('df1=%d,df2=%d,zeta.s=%g: nominal: %.2g%%; achieved: %.2g%%\n',df1,df2,zeta.s,
											#100*nominal.coverage,100*coverage))
					expect_equal(coverage,nominal.coverage,tolerance=alpha.tol,scale=1)
				}
			}
		}
	}
})#UNFOLD
#UNFOLD
context("estimation functions: prediction intervals")#FOLDUP
test_that("predint right output",{#FOLDUP
	set.char.seed("0919609e-4e99-42f2-b04c-89e2d2faaa4f")

	for (nasset in c(1,2,4,8)) {
		x <- matrix(rnorm(100*nasset),ncol=nasset)
		srx <- as.sr(x)
		aci <- predint(srx,oosdf=500)
		expect_equal(nrow(aci),nasset)
		expect_equal(ncol(aci),2)
		expect_true(is.matrix(aci))
	}

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("predint runs at all",{#FOLDUP
	set.char.seed("080c6f73-834e-4d10-a6fa-4b27dc266b24")

	ngen <- ceiling(THOROUGHNESS * 32)
	alpha.tol = 0.05 + 0.10 / THOROUGHNESS

	isn <- 200
	oosn <- 100

	ope <- 253
	sg <- 0.013
	for (nyr in c(3,6,9)) {
		df <- ceiling(ope * nyr)
		for (zeta in c(-1.0,0.0,1.0,2.0)) {
			x <- rnorm(isn,mean=(zeta/sqrt(ope))*sg,sd=sg)
			#y <- rnorm(oosn,mean=(zeta/sqrt(ope))*sg,sd=sg)
			for (nominal.coverage in c(0.90,0.95)) {
				aci <- predint(x,oosdf=oosn-1,ope=1,level=nominal.coverage)
				aci2 <- predint(as.sr(x),oosdf=oosn-1,ope=1,level=nominal.coverage)
				# no ope:
				noaci <- predint(x,oosdf=oosn-1,level=nominal.coverage)
			}
			# corner cases:
			iinf <- predint(x,oosdf=oosn-1,level.lo=0,level.hi=1)
			expect_true(all(is.infinite(iinf)))
		}
	}
	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("predint not fooled by annualization",{#FOLDUP
	set.char.seed("94cbb2a6-b497-48d4-9139-e987364f8a0f")

	isn <- 200
	oosn <- 100
	ope <- 253
	sg <- 0.013
	zeta <- 1
	x <- rnorm(isn,mean=(zeta/sqrt(ope))*sg,sd=sg)
	srx1 <- as.sr(x,ope=1)
	srx2 <- as.sr(x,ope=ope)
	nominal.coverage <- 0.90
	aci1 <- predint(srx1,oosdf=oosn-1,ope=1,level=nominal.coverage)
	aci2 <- predint(srx2,oosdf=oosn-1,ope=1,level=nominal.coverage)
	errs <- unlist(aci1) - unlist(aci2)

	expect_less_than(max(abs(errs)),1e-4)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
