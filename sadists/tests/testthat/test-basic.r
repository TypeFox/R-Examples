# /usr/bin/r
#
# Copyright 2015-2016 Steven E. Pav. All Rights Reserved.
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
#
# Created: 2016.01.04
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}




test_dpqr <- function(dpqr,nobs=1e4,...) {
	rv <- sort(dpqr$r(nobs,...))
	# stress test the d function...
	subrv <- rv[1:3]
	for (log.d in c(FALSE,TRUE)) {
		throwaway <- dpqr$d(subrv,...,log=log.d)
	}
	# stress test the p function...
	subrv <- rv[1:3]
	for (lower.tail in c(FALSE,TRUE)) { 
		for (log.p in c(FALSE,TRUE)) {
			throwaway <- dpqr$p(rv,...,lower.tail=lower.tail,log.p=log.p)
		}
	}
	# stress test the q function...
	somep <- ppoints(5)
	for (lower.tail in c(FALSE,TRUE)) { 
		for (log.p in c(FALSE,TRUE)) {
			if (log.p) { 
				checkp <- log(somep)
			} else {
				checkp <- somep
			}
			throwaway <- dpqr$q(checkp,...,lower.tail=lower.tail,log.p=log.p)
		}
	}
	
	# not sure what to do with these:
	#dv <- dpqr$d(rv,...)
	#pv <- dpqr$p(rv,...)
	#pvqv <- dpqr$q(pv,...)
	# dunno how to test the qv function for correctness.


	# anderson darling test: tests p function
	pp <- ppoints(nobs)
	qv <- dpqr$q(pp,...)
	S <- sum(((2 * seq(1,nobs) - 1) / nobs) * (log(pp) + log(1 - rev(pp))))
	A2 <- -nobs - S
	expect_lt(A2,1.0)

	# should check whether the p function and q function are near inverses ... 
	qvpv <- dpqr$p(qv,...)
	maxdev <- max(abs(sort(qvpv) - sort(pp)))
	expect_lt(maxdev,10.0/nobs)

}

context("Basic Operations")#FOLDUP

nobs_chk <- 2^11

test_that("sum of (non-central) chi-squares to a power",{
	set.char.seed('f8d988e6-e9ba-4448-aefe-51c50c0dbc02')
	wts <- c(-1,1,3,-3)
	df <- c(100,200,100,50)
	ncp <- c(0,1,0.5,2)
	pow <- c(1,0.5,2,1.5)
	test_dpqr(list(d=dsumchisqpow,p=psumchisqpow,q=qsumchisqpow,r=rsumchisqpow),nobs=nobs_chk,wts,df,ncp,pow)

	# sentinel:
	expect_true(TRUE)
})

test_that("K-prime distribution",{
	set.char.seed('0af8a714-9af6-4279-8422-b790755d8c89')
	v1 <- 50
	v2 <- 80
	a <- 0.5
	b <- 1.5
	test_dpqr(list(d=dkprime,p=pkprime,q=qkprime,r=rkprime),nobs=nobs_chk,v1,v2,a,b)

	set.char.seed('1d6669ae-0d9c-4662-9dbe-dfc705eee9c3')
	v1 <- 10
	v2 <- 30
	a <- 0
	b <- 1.5
	test_dpqr(list(d=dkprime,p=pkprime,q=qkprime,r=rkprime),nobs=nobs_chk,v1,v2,a,b)

	set.char.seed('11df12e2-8542-4cf0-adaf-250ff652bed2')
	v1 <- 100
	v2 <- 2
	a <- 2
	b <- 0
	test_dpqr(list(d=dkprime,p=pkprime,q=qkprime,r=rkprime),nobs=nobs_chk,v1,v2,a,b)

	# sentinel:
	expect_true(TRUE)
})

test_that("Lambda prime distribution",{
	set.char.seed('6799b96d-13d1-4733-afaa-bc7dd12e7c9e')

	df <- 50
	ts <- 1.5
	test_dpqr(list(d=dlambdap,p=plambdap,q=qlambdap,r=rlambdap),nobs=nobs_chk,df,ts)
	# sentinel
	expect_true(TRUE)
})

test_that("Upsilon distribution",{
	set.char.seed('ee481c35-006d-4be5-9989-1cecc3d8a2df')
	df <- c(30,50,100,20,10)
	ts <- c(-3,2,5,-4,1)
	test_dpqr(list(d=dupsilon,p=pupsilon,q=qupsilon,r=rupsilon),nobs=nobs_chk,df,ts)
	# sentinel
	expect_true(TRUE)
})
test_that("Doubly non-central F distribution",{
	set.char.seed('01116ab0-0961-44d8-9c02-3eea00856b4b')
	df1 <- 40
	df2 <- 80
	ncp1 <- 1.5
	ncp2 <- 2.5
	test_dpqr(list(d=ddnf,p=pdnf,q=qdnf,r=rdnf),nobs=nobs_chk,df1,df2,ncp1,ncp2)
	# sentinel
	expect_true(TRUE)
})
test_that("Doubly non-central t distribution",{
	set.char.seed('e3cf77e4-8fa1-4e79-b1fe-93307f4c097c')
	df <- 75
	ncp1 <- 2
	ncp2 <- 3
	test_dpqr(list(d=ddnt,p=pdnt,q=qdnt,r=rdnt),nobs=nobs_chk,df,ncp1,ncp2)
	# sentinel
	expect_true(TRUE)
})
test_that("Doubly non-central Beta distribution",{
	set.char.seed('b900579f-ae31-4e86-bdbe-7140e1454ea0')
	df1 <- 40
	df2 <- 80
	ncp1 <- 1.5
	ncp2 <- 2.5
	test_dpqr(list(d=ddnbeta,p=pdnbeta,q=qdnbeta,r=rdnbeta),nobs=nobs_chk,df1,df2,ncp1,ncp2)
	# sentinel
	expect_true(TRUE)
})
test_that("Doubly non-central Eta distribution",{
	set.char.seed('7bd4fe79-7323-4f2a-bbc0-beb44e9c5024')
	df <- 100
	ncp1 <- 0.5
	ncp2 <- 2.5
	test_dpqr(list(d=ddneta,p=pdneta,q=qdneta,r=rdneta),nobs=nobs_chk,df,ncp1,ncp2)
	# sentinel
	expect_true(TRUE)
})
test_that("Sum of logs of (non-central) chi-squares",{
	set.char.seed('fcc53464-e8c0-4b85-9404-65417bcb7740')
	wts <- c(5,-4,10,-15)
	df <- c(100,200,100,50)
	ncp <- c(0,1,0.5,2)
	test_dpqr(list(d=dsumlogchisq,p=psumlogchisq,q=qsumlogchisq,r=rsumlogchisq),nobs=nobs_chk,wts,df,ncp)
	# sentinel
	expect_true(TRUE)
})
test_that("Product of (non-central) chi-squares to a power",{
	set.char.seed('d269f717-4c62-4c1c-a8a9-4c67c5d5b8e8')
	df <- c(100,200,100,50)
	ncp <- c(0,1,0.5,2)
	pow <- c(1,0.5,2,1.5)
	test_dpqr(list(d=dprodchisqpow,p=pprodchisqpow,q=qprodchisqpow,r=rprodchisqpow),nobs=nobs_chk,df,ncp,pow)
	# sentinel
	expect_true(TRUE)
})
test_that("Product of doubly non-central F variates",{
	set.char.seed('7235cdf6-8b4c-4988-8034-d0fdba778bdf')
	df1 <- c(10,20,5)
	df2 <- c(1000,500,150)
	ncp1 <- c(1,0,2.5)
	ncp2 <- c(0,1.5,5)
	test_dpqr(list(d=dproddnf,p=pproddnf,q=qproddnf,r=rproddnf),nobs=nobs_chk,df1,df2,ncp1,ncp2)
	# sentinel
	expect_true(TRUE)
})

#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
