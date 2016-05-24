# /usr/bin/r
#
# Copyright 2015-2015 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of madness.
#
# madness is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# madness is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with madness.  If not, see <http://www.gnu.org/licenses/>.
#
# Created: 2015.11.22
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#library(expm)

set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}

# evaluate the error between a numerical approximation and a computed
# value, both of which may be erroneous?
errit <- function(apx,cmp,eps=1e-12) {
	merror <- abs(apx - cmp)
	rerror <- merror / pmax(eps^0.333,0.5 * (abs(apx) + abs(cmp)))
	rerror[(abs(apx) < eps^2) & (abs(cmp) < eps^2)] <- 0
	rerror
}

# now the harness
comp_err <- function(xval,thefun,scalfun=thefun,eps=1e-8) {
	xobj <- madness(val=xval,vtag='x',xtag='x')
	yobj <- thefun(xobj)
	# compute the error between the function applied to madness
	# and the function on the scalar value.
	ynum <- scalfun(xval)
	f_err <- errit(ynum,val(yobj),eps)

	# now the derivatives
	xval <- val(xobj)
  dapx <- numderiv(scalfun,xval,eps=eps,type='central')
	# compute error:
	dcmp <- dvdx(yobj)
	dim(dcmp) <- dim(dapx)

	d_err <- errit(dapx,dcmp,eps)
	
	retv <- max(max(abs(d_err)),max(abs(f_err)))
}

expect_small_err <- function(xval,thefun,scalfun=thefun,eps=1e-8,errtol=1e-6) {
	expect_less_than(comp_err(xval,thefun=thefun,scalfun=scalfun,eps=eps),errtol)
}

context("Basic Operations")#FOLDUP

test_that("arith functions",{#FOLDUP
	set.char.seed("dee9af9b-cb59-474f-ac3b-acd60faa8ba2")
	xval <- matrix(1 + runif(4*4),nrow=4)
	yval <- matrix(1 + runif(length(xval)),nrow=nrow(xval))
	expect_small_err(xval,function(x) { + x })
	expect_small_err(xval,function(x) { - x })

	expect_small_err(xval,function(x) { x + x })
	expect_small_err(xval,function(x) { x + yval })
	expect_small_err(xval,function(x) { yval + x })

	expect_small_err(xval,function(x) { x - x })
	expect_small_err(xval,function(x) { x - yval })
	expect_small_err(xval,function(x) { yval - x })

	expect_small_err(xval,function(x) { x * x })
	expect_small_err(xval,function(x) { x * yval })
	expect_small_err(xval,function(x) { yval * x })

	expect_small_err(xval,function(x) { x / x })
	expect_small_err(xval,function(x) { x / yval })
	expect_small_err(xval,function(x) { yval / x })

	expect_small_err(xval,function(x) { x ^ x })
	expect_small_err(xval,function(x) { x ^ yval })
	expect_small_err(xval,function(x) { yval ^ x })

	expect_small_err(xval,function(x) { x %*% x })
	expect_small_err(xval,function(x) { x %*% yval })
	expect_small_err(xval,function(x) { yval %*% x })

	expect_small_err(xval,function(x) { x %x% x })
	expect_small_err(xval,function(x) { x %x% yval })
	expect_small_err(xval,function(x) { yval %x% x })

	expect_small_err(xval,function(x) { crossprod(x) })
	expect_small_err(xval,function(x) { crossprod(x,x) })
	expect_small_err(xval,function(x) { crossprod(x,yval) })
	expect_small_err(xval,function(x) { crossprod(yval,x) })

	expect_small_err(xval,function(x) { tcrossprod(x) })
	expect_small_err(xval,function(x) { tcrossprod(x,x) })
	expect_small_err(xval,function(x) { tcrossprod(x,yval) })
	expect_small_err(xval,function(x) { tcrossprod(yval,x) })
	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("bind functions",{#FOLDUP
	set.char.seed("f459e4a4-2b1f-4902-9f5c-a78ee3302e96")
	xval <- matrix(1 + runif(4*4),nrow=4)
	yval <- matrix(1 + runif(length(xval)),nrow=nrow(xval))
	
	acol <- matrix(1 + runif(nrow(xval)),ncol=1)
	arow <- matrix(1 + runif(ncol(xval)),nrow=1)

	expect_small_err(xval,function(x) { cbind(x) })
	expect_small_err(xval,function(x) { cbind(x,x) })
	expect_small_err(xval,function(x) { cbind(x,x,x) })
	expect_small_err(xval,function(x) { cbind(x,x,x,x) })
	expect_small_err(xval,function(x) { cbind(x,yval) })
	expect_small_err(xval,function(x) { cbind(x,acol) })
	expect_small_err(xval,function(x) { cbind(x,acol,x) })
	expect_small_err(xval,function(x) { cbind(yval,x) })
	expect_small_err(xval,function(x) { cbind(x,yval,x) })
	expect_small_err(xval,function(x) { cbind(yval,x,yval) })

	expect_small_err(xval,function(x) { rbind(x) })
	expect_small_err(xval,function(x) { rbind(x,x) })
	expect_small_err(xval,function(x) { rbind(x,x,x) })
	expect_small_err(xval,function(x) { rbind(x,x,x,x) })
	expect_small_err(xval,function(x) { rbind(x,yval) })
	expect_small_err(xval,function(x) { rbind(x,arow) })
	expect_small_err(xval,function(x) { rbind(x,arow,x) })
	expect_small_err(xval,function(x) { rbind(yval,x) })
	expect_small_err(xval,function(x) { rbind(x,yval,x) })
	expect_small_err(xval,function(x) { rbind(yval,x,yval) })

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("elwise functions",{#FOLDUP
	set.char.seed("05ffaa40-902d-430a-a47f-63938b921306")
	xval <- matrix(1 + runif(4*4),nrow=4)

	expect_small_err(xval,function(x) { abs(x) })
	expect_small_err(xval,function(x) { exp(x) })
	expect_small_err(abs(xval),function(x) { log(x) })
	expect_small_err(abs(xval),function(x) { log10(x) })

	expect_small_err(xval,function(x) { sqrt(x) })
	expect_small_err(xval,function(x) { sin(x) })
	expect_small_err(xval,function(x) { cos(x) })
	expect_small_err(xval,function(x) { tan(x) },eps=1e-10)
	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("matwise functions",{#FOLDUP
	set.char.seed("1860a172-ba3e-4873-950a-75c9b771134d")
	zval <- matrix(0.01 + runif(4*100,min=0,max=0.05),nrow=4)
	xval <- tcrossprod(zval)

	expect_small_err(xval,function(x) { sqrtm(x) })
	#expect_small_err(xval,function(x) { logm(x) })
	#expect_small_err(xval,function(x) { expm(x) })

	# the 'wrong way to test'
	expect_gt(comp_err(xval,function(x) { chol(x) }),1e-6)

	# chol has hidden symmetry:
	fsym <- function(x) { 0.5 * (x + t(x)) }
	expect_small_err(xval,function(x) { chol(fsym(x)) })
	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("sum functions",{#FOLDUP
	set.char.seed("25e43832-3030-40cf-acf8-fa43bd56dc09")
	xval <- matrix(1 + runif(4*4),nrow=4)

	expect_small_err(xval,function(x) { matrix.trace(x) },errtol=1e-6)

	for (na.rm in c(FALSE,TRUE)) {
		expect_small_err(xval,function(x) { colSums(x,na.rm=na.rm) })
		expect_small_err(xval,function(x) { colMeans(x,na.rm=na.rm) })
		expect_small_err(xval,function(x) { rowSums(x,na.rm=na.rm) })
		expect_small_err(xval,function(x) { rowMeans(x,na.rm=na.rm) })
	}

	xval <- matrix(runif(6*6,min=1,max=2),nrow=6)
	xval[xval < 1.2] <- NA
	na.rm <- TRUE
	expect_small_err(xval,function(x) { colSums(x,na.rm=na.rm) })
	expect_small_err(xval,function(x) { colMeans(x,na.rm=na.rm) })
	expect_small_err(xval,function(x) { rowSums(x,na.rm=na.rm) })
	expect_small_err(xval,function(x) { rowMeans(x,na.rm=na.rm) })
	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("sumprod functions",{#FOLDUP
	set.char.seed("a5392bf9-9b9b-4411-87f5-9e40365c73d9")
	xval <- matrix(runif(4*4,min=1,max=2),nrow=4)
	for (na.rm in c(FALSE,TRUE)) {
		expect_small_err(xval,function(x) { sum(x,na.rm=na.rm) })
		expect_small_err(xval,function(x) { prod(x,na.rm=na.rm) })
	}

	xval <- matrix(runif(6*6,min=1,max=2),nrow=6)
	xval[xval < 1.2] <- NA
	na.rm <- TRUE
	expect_small_err(xval,function(x) { sum(x,na.rm=na.rm) })
	expect_small_err(xval,function(x) { prod(x,na.rm=na.rm) })

	xval <- matrix(rnorm(6*6),nrow=6)
	xval[xval < -1] <- NA
	na.rm <- TRUE
	expect_small_err(xval,function(x) { sum(x,na.rm=na.rm) })
	expect_small_err(xval,function(x) { prod(x,na.rm=na.rm) })
	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("maxmin functions",{#FOLDUP
	set.char.seed("c52e167b-14ea-4827-ad47-5ee165929ed7")
	xval <- matrix(runif(4*4,min=1,max=2),nrow=4)
	for (na.rm in c(FALSE,TRUE)) {
		expect_small_err(xval,function(x) { max(x,na.rm=na.rm) })
		expect_small_err(xval,function(x) { min(x,na.rm=na.rm) })
	}

	xval <- matrix(runif(6*6,min=1,max=2),nrow=6)
	xval[xval < 1.2] <- NA
	na.rm <- TRUE
	expect_small_err(xval,function(x) { max(x,na.rm=na.rm) })
	expect_small_err(xval,function(x) { min(x,na.rm=na.rm) })

	xval <- matrix(rnorm(6*6),nrow=6)
	xval[xval < -1] <- NA
	na.rm <- TRUE
	expect_small_err(xval,function(x) { max(x,na.rm=na.rm) })
	expect_small_err(xval,function(x) { min(x,na.rm=na.rm) })
	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("outer functions",{#FOLDUP
	set.char.seed("354455b9-2b3f-40fb-b8e7-21a1302b48de")
	xval <- array(1+runif(2*3*4),dim=c(2,3,4))
	yval <- array(1+runif(3*2),dim=c(3,2))

	expect_small_err(xval,function(x) { outer(x,x,FUN='*') })
	expect_small_err(xval,function(x) { outer(x,x,FUN='+') })
	expect_small_err(xval,function(x) { outer(x,x,FUN='-') })
	expect_small_err(xval,function(x) { outer(x,x,FUN='/') })

	expect_small_err(xval,function(x) { outer(x,yval,FUN='*') })
	expect_small_err(xval,function(x) { outer(x,yval,FUN='+') })
	expect_small_err(xval,function(x) { outer(x,yval,FUN='-') })
	expect_small_err(xval,function(x) { outer(x,yval,FUN='/') })

	expect_small_err(xval,function(x) { outer(yval,x,FUN='*') })
	expect_small_err(xval,function(x) { outer(yval,x,FUN='+') })
	expect_small_err(xval,function(x) { outer(yval,x,FUN='-') })
	expect_small_err(xval,function(x) { outer(yval,x,FUN='/') })
	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("determinants",{#FOLDUP
	set.char.seed("081849a0-ab28-42ac-8d18-1963cb8a9a0a")
	xval <- matrix(1 + runif(4*4),nrow=4)

	# fuck det.
	#expect_small_err(xval,function(x) { det(x) })
	expect_small_err(xval,function(x) { dt <- determinant(x,logarithm=FALSE) ; dt$modulus },errtol=1e-6)
	expect_small_err(xval,function(x) { dt <- determinant(x,logarithm=TRUE)  ; dt$modulus },errtol=1e-6)
	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("reshape functions",{#FOLDUP
	set.char.seed("3d06aea4-c339-4630-a8db-3d56d6b6b687")
	xval <- matrix(1 + runif(4*4),nrow=4)
	xvec <- array(1 + runif(4*4),dim=c(16,1))

	expect_small_err(xval,function(x) { t(x) },errtol=1e-6)

	expect_small_err(xval,function(x) { vec(x) },errtol=1e-6)
	expect_small_err(xval,function(x) { vech(x) },errtol=1e-6)
	expect_small_err(xval,function(x) { vech(x,1) },errtol=1e-6)
	expect_small_err(xval,function(x) { vech(x,-1) },errtol=1e-6)

	expect_small_err(xval,function(x) { diag(x) },errtol=1e-6)
	
	expect_small_err(xval,function(x) { dim(x) <- c(prod(dim(x)),1); x },errtol=1e-6)
	expect_small_err(xval,function(x) { x[1,1,drop=FALSE] },errtol=1e-6)
	expect_small_err(xvec,function(x) { dim(x) <- c(prod(dim(x)),1); x },errtol=1e-6)
	expect_small_err(xvec,function(x) { x[1,1,drop=FALSE] },errtol=1e-6)

	expect_small_err(xval,function(x) { todiag(x) }, function(x) { diag(as.numeric(x)) },errtol=1e-6)

	expect_small_err(xval,function(x) { tril(x) }, function(x) { x[row(x) < col(x)] <- 0; x },errtol=1e-6)
	expect_small_err(xval,function(x) { triu(x) }, function(x) { x[row(x) > col(x)] <- 0; x },errtol=1e-6)
	

	xval <- array(1 + runif(2*3*4*5),dim=c(2,3,4,5))
	expect_small_err(xval,function(x) { aperm(x) },errtol=1e-6)
	expect_small_err(xval,function(x) { aperm(x,c(2,1,3,4)) },errtol=1e-6)
	expect_small_err(xval,function(x) { aperm(x,c(4,3,2,1)) },errtol=1e-6)

	# need better tests of these!
	xval <- array(1 + runif(2*3*4*5),dim=c(2,3,4,5))
	xdim <- dim(xval)
	expect_small_err(xval,function(x) { blockrep(x,c(1)) },function(x) { x },errtol=1e-6)
	expect_small_err(xval,function(x) { repto(x,xdim) },function(x) { x },errtol=1e-6)

	xval <- array(1 + runif(3*7),dim=c(3,7))
	xdim <- dim(xval)
	expect_small_err(xval,function(x) { blockrep(x,c(2,1)) },function(x) { rbind(x,x) },errtol=1e-6)
	expect_small_err(xval,function(x) { blockrep(x,c(3,1)) },function(x) { rbind(x,x,x) },errtol=1e-6)
	expect_small_err(xval,function(x) { blockrep(x,c(5,1)) },function(x) { rbind(x,x,x,x,x) },errtol=1e-6)
	expect_small_err(xval,function(x) { blockrep(x,c(1,2)) },function(x) { cbind(x,x) },errtol=1e-6)
	expect_small_err(xval,function(x) { blockrep(x,c(1,4)) },function(x) { cbind(x,x,x,x) },errtol=1e-6)
	expect_small_err(xval,function(x) { blockrep(x,c(2,3)) },function(x) { x2 <- rbind(x,x) ; cbind(x2,x2,x2) },errtol=1e-6)
	expect_small_err(xval,function(x) { blockrep(x,c(5,2)) },function(x) { x5 <- rbind(x,x,x,x,x) ; cbind(x5,x5) },errtol=1e-6)
	# now alternate dimensions?

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("vech/ivech functions",{#FOLDUP
	set.char.seed("a61a011f-d13c-4a70-8ab7-7003edbf65a0")
	xval <- matrix(1 + runif(4*4),nrow=4)
	xvec <- array(1 + runif(4*4),dim=c(16,1))

	expect_small_err(xval,function(x) { vec(x) },errtol=1e-6)
	expect_small_err(xval,function(x) { vech(x) },errtol=1e-6)
	expect_small_err(xval,function(x) { vech(x,1) },errtol=1e-6)
	expect_small_err(xval,function(x) { vech(x,-1) },errtol=1e-6)

	set.char.seed("815ac0d5-3e63-45cb-bf65-cb2bb26bb8d6")

	# ivech
	for (nr in c(3,6,10)) {
		xval <- matrix(1 + runif(nr),nrow=nr)
		for (symmetric in c(FALSE,TRUE)) {
			expect_small_err(xval,function(x) { ivech(x,0,symmetric=symmetric) },errtol=1e-6)
			expect_small_err(xval,function(x) { ivech(x,-1,symmetric=symmetric) },errtol=1e-6)
			expect_small_err(xval,function(x) { ivech(x,-2,symmetric=symmetric) },errtol=1e-6)
		}
	}
	for (nr in c(8,13)) {
		xval <- matrix(1 + runif(nr),nrow=nr)
		expect_small_err(xval,function(x) { ivech(x,1,symmetric=FALSE) },errtol=1e-6)
		expect_small_err(xval,function(x) { ivech(x,1,symmetric=FALSE) },errtol=1e-6)
	}
	for (nr in c(15,22)) {
		xval <- matrix(1 + runif(nr),nrow=nr)
		expect_small_err(xval,function(x) { ivech(x,2,symmetric=FALSE) },errtol=1e-6)
		expect_small_err(xval,function(x) { ivech(x,2,symmetric=FALSE) },errtol=1e-6)
	}
	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("solve functions",{#FOLDUP
	set.char.seed("232ba1a1-7751-40be-866f-a8e2122c2ace")
	prex <- matrix(1 + runif(100*8),ncol=8)
	xval <- crossprod(prex)
	yval <- array(1 + runif(nrow(xval)),dim=c(nrow(xval),1))

	expect_small_err(xval,function(x) { solve(x) },eps=1e-6,errtol=1e-6)
	expect_small_err(xval,function(x) { solve(x,yval) },eps=1e-6,errtol=1e-6)
	expect_small_err(xval,function(x) { solve(x,x[,1,drop=FALSE]) },eps=1e-6,errtol=1e-5)
	expect_small_err(xval,function(x) { solve(xval,x[,1,drop=FALSE]) },eps=1e-6,errtol=1e-5)
	expect_small_err(xval,function(x) { solve(x,as.numeric(yval)) },eps=1e-6,errtol=1e-6)

	# numeric left only works in 1d case. so be degenerate
	xval <- array(rnorm(1),dim=c(1,1))
	yval <- runif(1)
	expect_small_err(xval,function(x) { solve(yval,x) },eps=1e-6,errtol=1e-6)

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("eigen functions",{#FOLDUP
	set.char.seed("c0529b18-8873-4353-9c41-b4c8bc24d319")
	xval <- matrix(1 + runif(4*4),nrow=4)
	xval <- xval + t(xval)

	expect_small_err(xval,function(x) { eigen(0.5 * (x + t(x)),symmetric=TRUE)$values },eps=1e-7,errtol=1e-6)
	expect_small_err(xval,function(x) { 
										 ev <- eigen(0.5 * (x + t(x)),symmetric=TRUE)
										 # force first row to be all positive
										 ev$vectors %*% diag(sign(val(ev$vectors)[1,,drop=TRUE]))
									 },
									 scalfun = function(x) {
										 ev <- eigen(0.5 * (x + t(x)),symmetric=TRUE)
										 # force first row to be all positive
										 ev$vectors %*% diag(sign(ev$vectors[1,,drop=TRUE]))
									 }, eps=1e-7,errtol=1e-6)

	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("norm functions",{#FOLDUP
	set.char.seed("e74da7ce-92f1-41ee-96cd-fe8201da753f")
	xval <- matrix(1 + runif(4*4),nrow=4)

	expect_small_err(xval,function(x) { norm(x) },errtol=1e-6)

	expect_small_err(xval,function(x) { norm(x,'O') },errtol=1e-6)
	expect_small_err(xval,function(x) { norm(x,'o') },errtol=1e-6)
	expect_small_err(xval,function(x) { norm(x,'1') },errtol=1e-6)
	expect_small_err(xval,function(x) { norm(x,'I') },errtol=1e-6)
	expect_small_err(xval,function(x) { norm(x,'i') },errtol=1e-6)
	expect_small_err(xval,function(x) { norm(x,'M') },errtol=1e-6)
	expect_small_err(xval,function(x) { norm(x,'m') },errtol=1e-6)
	expect_small_err(xval,function(x) { norm(x,'F') },errtol=1e-6)
	expect_small_err(xval,function(x) { norm(x,'f') },errtol=1e-6)
# Matrix::norm does not support type '2'
	expect_small_err(xval,function(x) { norm(x,'2') },function(x) { base::norm(x,'2') },errtol=1e-6)

	expect_small_err(xval,function(x) { maxeig(x) },
																function(x) { 
																	usv <- svd(x,1,1)
																	as.numeric(usv$d[1])
																},errtol=1e-6)
	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD

#UNFOLD

context("Requiring Chain Rule")#FOLDUP

test_that("round one",{#FOLDUP
	set.char.seed("dee9af9b-cb59-474f-ac3b-acd60faa8ba2")
	xval <- matrix(1 + runif(4*4),nrow=4)
	yval <- matrix(1 + runif(length(xval)),nrow=nrow(xval))
	expect_small_err(xval,function(x) { norm(crossprod(x),'O') },eps=1e-7,errtol=1e-5)
	expect_small_err(xval,function(x) { norm(crossprod(x^x),'M') },eps=1e-7,errtol=1e-5)
	expect_small_err(xval,function(x) { norm(abs(x) %*% t(x),'I') },eps=1e-7,errtol=1e-5)

	expect_small_err(xval,function(x) { tcrossprod(sin(x)) },eps=1e-07,errtol=1e-5)
	expect_small_err(xval,function(x) { cos(crossprod(sin(x))) },eps=1e-07,errtol=1e-5)

	# sentinel:
	expect_true(TRUE)
})#UNFOLD

#UNFOLD

context("Inverses")#FOLDUP

test_that("vech",{#FOLDUP
	set.char.seed("dee9af9b-cb59-474f-ac3b-acd60faa8ba2")
	xval <- matrix(1 + runif(4*4),nrow=4)
	yval <- vech(xval,0)
	xival <- ivech(yval,0)
	yival <- vech(xival,0)

	expect_less_than(max(errit(yval,yival)),1e-12)

	xsym <- xval + t(xval)
	yval <- vech(xsym,0)
	xsval <- ivech(yval,0,symmetric=TRUE)
	expect_less_than(max(errit(xsym,xsval)),1e-12)

	# sentinel:
	expect_true(TRUE)
})#UNFOLD

#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
