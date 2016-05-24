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
# Created: 2013.08.14
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
test_that("sr_test",{#FOLDUP
	set.char.seed("0b144107-4de8-4e00-95f7-d746db3aef8e")
	X <- matrix(rnorm(1000*1),ncol=1)
	Y <- matrix(rnorm(length(X)),ncol=1)
	Y2 <- matrix(rnorm(2*length(X)),ncol=1)

	fooz <- sr_test(X,alternative="two.sided")
	fooz <- sr_test(X,alternative="less")
	fooz <- sr_test(X,alternative="greater")

	fooz <- sr_test(as.sr(X),alternative="two.sided")

	# with a zeta
	zeta <- 1.0
	fooz <- sr_test(X,zeta=zeta,alternative="two.sided")
	fooz <- sr_test(X,zeta=zeta,alternative="less")
	fooz <- sr_test(X,zeta=zeta,alternative="greater")

	# and an ope
	ope <- 252
	fooz <- sr_test(X,zeta=zeta,ope=ope,alternative="two.sided")
	fooz <- sr_test(X,zeta=zeta,ope=ope,alternative="less")
	fooz <- sr_test(X,zeta=zeta,ope=ope,alternative="greater")

	# X and Y
	fooz <- sr_test(X,Y,ope=ope,paired=TRUE,alternative="two.sided")
	fooz <- sr_test(X,Y,ope=ope,paired=TRUE,alternative="less")
	fooz <- sr_test(X,Y,ope=ope,paired=TRUE,alternative="greater")

	fooz <- sr_test(X,Y,ope=ope,paired=FALSE,alternative="two.sided")
	fooz <- sr_test(X,Y,ope=ope,paired=FALSE,alternative="less")
	fooz <- sr_test(X,Y,ope=ope,paired=FALSE,alternative="greater")

	fooz <- sr_test(X,Y2,ope=ope,paired=FALSE,alternative="two.sided")
	fooz <- sr_test(X,Y2,ope=ope,paired=FALSE,alternative="less")
	fooz <- sr_test(X,Y2,ope=ope,paired=FALSE,alternative="greater")

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("sr_unpaired_test",{#FOLDUP
	set.char.seed("6036d6e3-b3e9-415a-bd2b-bfc19dbd45ea")
	X <- matrix(rnorm(1000*6),ncol=6)
	inp <- as.sr(X)
	nv <- dim(X)[2]
	etc <- sr_unpaired_test(inp)
	etc <- sr_unpaired_test(inp,alternative='two.sided')
	etc <- sr_unpaired_test(inp,alternative='less')
	etc <- sr_unpaired_test(inp,alternative='greater')
	etc <- sr_unpaired_test(inp,ope=1)
	etc <- sr_unpaired_test(inp,ope=10)
	etc <- sr_unpaired_test(inp,contrasts=runif(nv))
	etc <- sr_unpaired_test(inp,contrasts=runif(nv),null.value=0.5)
	etc <- sr_unpaired_test(inp,contrasts=runif(nv),null.value=0.5,alternative='less')
	etc <- sr_unpaired_test(inp,contrasts=runif(nv),null.value=0.5,alternative='greater')

	inp <- apply(X,2,as.sr)
	nv <- length(inp)
	etc <- sr_unpaired_test(inp)
	etc <- sr_unpaired_test(inp,alternative='two.sided')
	etc <- sr_unpaired_test(inp,alternative='less')
	etc <- sr_unpaired_test(inp,alternative='greater')
	etc <- sr_unpaired_test(inp,ope=1)
	etc <- sr_unpaired_test(inp,ope=10)
	etc <- sr_unpaired_test(inp,contrasts=runif(nv))
	etc <- sr_unpaired_test(inp,contrasts=runif(nv),null.value=0.5)
	etc <- sr_unpaired_test(inp,contrasts=runif(nv),null.value=0.5,alternative='less')
	etc <- sr_unpaired_test(inp,contrasts=runif(nv),null.value=0.5,alternative='greater')

	# bad units, complain
	X1 <- rnorm(100)
	X2 <- rnorm(200)
	inp1 <- as.sr(X1,ope=52)
	inp2 <- as.sr(X2,ope=12)
	expect_warning(sr_unpaired_test(list(inp1,inp2)))

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("sr_equality_test",{#FOLDUP
	set.char.seed("0b144107-4de8-4e00-95f7-d746db3aef8e")
	X <- matrix(rnorm(1000*5),ncol=5)
	Con = matrix(rnorm(dim(X)[2]),nrow=1)
	fooz <- sr_equality_test(X,type="chisq")
	fooz <- sr_equality_test(X,type="F")
	fooz <- sr_equality_test(X,type="t",contrasts=Con)
	fooz <- sr_equality_test(X,type="t",contrasts=Con,alternative='less')
	fooz <- sr_equality_test(X,type="t",contrasts=Con,alternative='greater')

	expect_warning(sr_equality_test(X,type="F",alternative='less'))
	expect_warning(sr_equality_test(X,type="F",alternative='greater'))

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("sropt_test",{#FOLDUP
	set.char.seed("02a77746-3f08-4320-9696-46521ab4de37")
	X <- matrix(rnorm(1000*10),ncol=10)

	fooz <- sropt_test(X,alternative="two.sided")
	fooz <- sropt_test(X,alternative="less")
	fooz <- sropt_test(X,alternative="greater")

	fooz <- sropt_test(as.sropt(X),alternative="two.sided")
	fooz <- sropt_test(as.sropt(X),alternative="less")
	fooz <- sropt_test(as.sropt(X),alternative="greater")

	# with a zeta
	zeta.s <- 1.0
	fooz <- sropt_test(X,zeta.s=zeta.s,alternative="two.sided")
	fooz <- sropt_test(X,zeta.s=zeta.s,alternative="less")
	fooz <- sropt_test(X,zeta.s=zeta.s,alternative="greater")

	# and an ope
	ope <- 252
	fooz <- sropt_test(X,zeta.s=zeta.s,ope=ope,alternative="two.sided")
	fooz <- sropt_test(X,zeta.s=zeta.s,ope=ope,alternative="less")
	fooz <- sropt_test(X,zeta.s=zeta.s,ope=ope,alternative="greater")

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("power_sr_test",{#FOLDUP
	set.char.seed("5e24476e-4001-472d-bd63-8a660c3737b4")

	ope <- 252
	apo <- power.sr_test(n=1000,zeta=1.0,sig.level=0.05,alternative="one.sided",ope=ope)
	apo <- power.sr_test(n=1000,power=0.3,sig.level=0.05,alternative="one.sided",ope=ope)
	# this is dumb:
	#apo <- power.sr_test(n=1000,zeta=0.9,power=0.3,alternative="one.sided",ope=ope)
	apo <- power.sr_test(n=1000,zeta=0.9,power=0.3,sig.level=NULL,alternative="one.sided",ope=ope)

	apo <- power.sr_test(n=1000,zeta=1.0,sig.level=0.05,alternative="two.sided",ope=ope)
	apo <- power.sr_test(n=1000,power=0.3,sig.level=0.05,alternative="two.sided",ope=ope)
	# this is dumb:
	#apo <- power.sr_test(n=1000,zeta=0.9,power=0.3,alternative="two.sided",ope=ope)
	apo <- power.sr_test(n=1000,zeta=0.9,power=0.3,sig.level=NULL,alternative="two.sided",ope=ope)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("power_sropt_test",{#FOLDUP
	set.char.seed("d486234e-e897-4dc9-ae8b-c4312b48c8b4")

	ope <- 1
	apo <- power.sropt_test(df1=10,df2=1000,zeta.s=1.0,sig.level=0.05,power=NULL,ope=ope)
	apo <- power.sropt_test(df1=10,df2=1000,zeta.s=2.0,sig.level=NULL,power=0.9,ope=ope)
	apo <- power.sropt_test(df1=10,df2=1000,zeta.s=NULL,sig.level=0.05,power=0.9,ope=ope)
	#apo <- power.sropt_test(df1=10,df2=NULL,zeta.s=1.0,sig.level=0.05,power=0.3,ope=ope)
	#apo <- power.sropt_test(df1=NULL,df2=1000,zeta.s=1.0,sig.level=0.05,power=0.4,ope=ope)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD
context("test hypothesis tests")#FOLDUP
test_that("sr_equality_test uniformity",{#FOLDUP
	ngen <- ceiling(THOROUGHNESS * 32)
	alpha.floor = 0.001 + 0.003 * (THOROUGHNESS / (1 + THOROUGHNESS))

	set.char.seed("0b144107-4de8-4e00-95f7-d746db3aef8e")
	ope <- 253
	for (nyr in c(2,4)) {
		nday <- ceiling(ope * nyr)
		for (nstock in c(2,4)) {
			for (mu in c(0,0.1)) {
				pvs <- rep(0,ngen)
				for (iii in (1:ngen)) {
					X <- matrix(rnorm(nday * nstock,mean=mu),ncol=nstock)
					etc <- sr_equality_test(X)
					pvs[iii] <- etc$p.value
				}
				meta.test <- ks.test(pvs,punif)
				expect_true(meta.test$p.value > alpha.floor)
			}
		}
	}

	if (require(sandwich)) {
		ngen <- ceiling(THOROUGHNESS * 16)

		set.char.seed("d1cbe090-dca6-4543-8332-1da96f1665a8")
		ope <- 253
		for (nyr in c(1,2)) {
			nday <- ceiling(ope * nyr)
			for (nstock in c(2)) {
				for (mu in c(0,0.1)) {
					pvs <- rep(0,ngen)
					for (iii in (1:ngen)) {
						X <- matrix(rnorm(nday * nstock,mean=mu),ncol=nstock)
						etc <- sr_equality_test(X,vcov.func=vcovHAC)
						pvs[iii] <- etc$p.value
					}
					meta.test <- ks.test(pvs,punif)
					expect_true(meta.test$p.value > alpha.floor)
				}
			}
		}
	}


})#UNFOLD
test_that("sr_test one sample uniformity under null",{#FOLDUP
	ngen <- ceiling(THOROUGHNESS * 64)
	alpha.floor = 0.001 + 0.003 * (THOROUGHNESS / (1 + THOROUGHNESS))

	set.char.seed("3e4249f2-18d4-4d5c-98be-c32cc05e74c7")
	ope <- 253
	for (nyr in c(2,4)) {
		nday <- ceiling(ope * nyr)
		for (psi in c(0,0.1)) {
			for (sg in c(0.01,0.02)) {
				pvs <- rep(0,ngen)
				for (iii in (1:ngen)) {
					X <- rnorm(nday,mean=psi*sg,sd=sg)
					etc <- sr_test(X,zeta=psi)
					pvs[iii] <- etc$p.value
				}
				meta.test <- ks.test(pvs,punif)
				expect_true(meta.test$p.value > alpha.floor)
			}
		}
	}
})#UNFOLD
test_that("sr_test two sample unpaired uniformity under null",{#FOLDUP
	ngen <- ceiling(THOROUGHNESS * 64)
	alpha.floor = 0.001 + 0.003 * (THOROUGHNESS / (1 + THOROUGHNESS))

	set.char.seed("43af0f33-ef6b-4c1a-ba78-a3e4ec190644")
	ope <- 253
	for (nyr in c(2,4)) {
		nday <- ceiling(ope * nyr)
		for (psi in c(0,0.1)) {
			for (sg in c(0.01,0.02)) {
				pvs <- rep(0,ngen)
				for (iii in (1:ngen)) {
					x <- rnorm(nday,mean=psi*sg,sd=sg)
					y <- rnorm(nday+20,mean=psi*sg,sd=sg)
					etc <- sr_test(x,y,paired=FALSE)
					pvs[iii] <- etc$p.value
				}
				meta.test <- ks.test(pvs,punif)
				expect_true(meta.test$p.value > alpha.floor)
			}
		}
	}
})#UNFOLD
test_that("sr_test two sample paired uniformity under null",{#FOLDUP
	ngen <- ceiling(THOROUGHNESS * 64)
	alpha.floor = 0.001 + 0.003 * (THOROUGHNESS / (1 + THOROUGHNESS))

	set.char.seed("68d2825b-e645-4bbe-94d0-566cd7d4d8e1")
	ope <- 253
	for (nyr in c(2,4)) {
		nday <- ceiling(ope * nyr)
		for (psi in c(0,0.1)) {
			for (sg in c(0.01,0.02)) {
				pvs <- rep(0,ngen)
				for (iii in (1:ngen)) {
					x <- rnorm(nday,mean=psi*sg,sd=sg)
					y <- rnorm(nday,mean=psi*sg,sd=sg)
					etc <- sr_test(x,y,paired=TRUE)
					pvs[iii] <- etc$p.value
				}
				meta.test <- ks.test(pvs,punif)
				expect_true(meta.test$p.value > alpha.floor)
			}
		}
	}
})#UNFOLD

test_that("sr_test one sample power under alternative",{#FOLDUP
	ngen <- ceiling(THOROUGHNESS * 32)
	alpha.floor = 0.001 + 0.003 * (THOROUGHNESS / (1 + THOROUGHNESS))

	set.char.seed("56a7e85d-516b-47b4-a98f-08127b57d13d")
	ope <- 253
	for (nyr in c(2,4)) {
		nday <- ceiling(ope * nyr)
		for (psi in c(0,0.1)) {
			for (sg in c(0.01,0.02)) {
				pvs.lo <- rep(0,ngen)
				pvs.hi <- rep(0,ngen)
				for (iii in (1:ngen)) {
					X <- rnorm(nday,mean=psi*sg,sd=sg)
					etc <- sr_test(X,zeta=psi-0.2,alternative='greater')
					pvs.lo[iii] <- etc$p.value
					etc <- sr_test(X,zeta=psi-0.2,alternative='less')
					pvs.hi[iii] <- etc$p.value
				}
				# these should imply each other, but it is a test.
				expect_true(max(pvs.lo) < 0.5)
				expect_true(min(pvs.hi) > 0.5)
			}
		}
	}
})#UNFOLD
test_that("sr_test two sample unpaired power under alternative",{#FOLDUP
	ngen <- ceiling(THOROUGHNESS * 32)
	alpha.floor = 0.001 + 0.003 * (THOROUGHNESS / (1 + THOROUGHNESS))

	set.char.seed("1ff5561a-90d7-4f95-a0a6-a60c8aaa31f5")
	ope <- 253
	for (nyr in c(2,4)) {
		nday <- ceiling(ope * nyr)
		for (psi in c(0,0.1)) {
			for (sg in c(0.01,0.02)) {
				pvs.lo <- rep(0,ngen)
				pvs.hi <- rep(0,ngen)
				for (iii in (1:ngen)) {
					x <- rnorm(nday,mean=psi*sg,sd=sg)
					y <- rnorm(nday+20,mean=(psi-0.2)*sg,sd=sg)
					etc <- sr_test(x,y,paired=FALSE,alternative='greater')
					pvs.lo[iii] <- etc$p.value
					etc <- sr_test(x,y,paired=FALSE,alternative='less')
					pvs.hi[iii] <- etc$p.value
				}
				# these should imply each other, but it is a test.
				expect_true(max(pvs.lo) < 0.5)
				expect_true(min(pvs.hi) > 0.5)
			}
		}
	}
})#UNFOLD
test_that("sr_test two sample paired power under alternative",{#FOLDUP
	ngen <- ceiling(THOROUGHNESS * 64)
	alpha.floor = 0.001 + 0.003 * (THOROUGHNESS / (1 + THOROUGHNESS))

	set.char.seed("b3e4d157-a618-4a8d-8556-600566afc741")
	ope <- 253
	for (nyr in c(2,4)) {
		nday <- ceiling(ope * nyr)
		for (psi in c(0,0.1)) {
			for (sg in c(0.01,0.02)) {
				pvs.lo <- rep(0,ngen)
				pvs.hi <- rep(0,ngen)
				for (iii in (1:ngen)) {
					x <- rnorm(nday,mean=psi*sg,sd=sg)
					y <- rnorm(nday,mean=(psi-0.25)*sg,sd=sg)
					etc <- sr_test(x,y,paired=TRUE,alternative='greater')
					pvs.lo[iii] <- etc$p.value
					etc <- sr_test(x,y,paired=TRUE,alternative='less')
					pvs.hi[iii] <- etc$p.value
				}
				# these should imply each other, but it is a test.
				expect_true(max(pvs.lo) < 0.5)
				expect_true(min(pvs.hi) > 0.5)
				expect_true(max(abs(pvs.lo + pvs.hi - 1.0)) < 0.01)
			}
		}
	}
})#UNFOLD

#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
