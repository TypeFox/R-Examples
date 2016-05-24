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
# Created: 2015.12.02
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

library(expm)

set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}

context("Basic Operations")#FOLDUP

test_that("initialize",{#FOLDUP
	set.char.seed("dcb575af-c16c-4bc9-a57e-f8f67b45462d")
	yt <- 'any'
	xt <- 'anx'

	xy <- data.frame(x=rnorm(100),y=runif(100),z=runif(100))
	amod <- lm(z ~ x + y,xy)
	amad <- as.madness(amod,xtag=xt,vtag=yt)
	amad <- as.madness(amod,xtag=xt)
	amad <- as.madness(amod,vtag=yt)
	amad <- as.madness(amod)

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("initialize errors",{#FOLDUP
	set.char.seed("ea2ca251-7b91-4835-a030-610f1835e995")

	yt <- 'any'
	xt <- 'anx'

	xval <- matrix(1 + runif(4*4),nrow=4)
	# make this the wrong size:
	ddd <- matrix(rnorm((1+length(xval))*5),ncol=5)
	# make this the wrong size
	vvv <- crossprod(matrix(rnorm(100*(1+ncol(ddd))),ncol=(1+ncol(ddd))))

	expect_error(dumb <- madness(xval,ddd,vtag=yt,xtag=xt,varx=vvv))
	expect_error(dumb <- madness(xval,ddd,vtag=yt,xtag=xt))
	expect_error(dumb <- madness(xval,diag(length(xval)),vtag=yt,xtag=xt,varx=vvv))
	expect_error(dumb <- madness(xval,'bogus',vtag=yt,xtag=xt))
	expect_error(dumb <- madness(xval,diag(length(xval)),vtag=23))
	expect_error(dumb <- madness(xval,diag(length(xval)),vtag=yt,xtag=17))
	expect_error(dumb <- madness(xval,diag(length(xval)),vtag=yt,xtag=xt,varx='dumby'))
# make this error out:
	#expect_error(dumb <- madness('not a numeric'))

	vvv <- matrix(rnorm(5*4),ncol=5)
	expect_error(dumb <- madness(xval,diag(length(xval)),vtag=yt,xtag=xt,varx=vvv))

# these are warnings
	xval <- rnorm(5)
	expect_warning(dumb <- madness(xval,vtag=yt,xtag=xt))
# fix this: should be a warning...
	xval <- matrix(1 + runif(4*4),nrow=4)
	nover <- 2
	ddd <- rnorm(nover*length(xval))
	expect_warning(dumb <- madness(xval,ddd,vtag=yt,xtag=xt))
	expect_warning(dumb <- madness(xval,ddd))
	expect_warning(dumb <- madness(xval,ddd,vtag=yt))
	expect_warning(dumb <- madness(xval,ddd,varx=diag(nover)))
	expect_error(dumb <- madness(xval,ddd,varx=diag(nover+1)))


	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("basic getters and setters",{#FOLDUP
	set.char.seed("dee9af9b-cb59-474f-ac3b-acd60faa8ba2")
	xval <- matrix(1 + runif(4*4),nrow=4)
	xmad <- madness(xval)

	expect_equal(xval,val(xmad))
	expect_equal(dim(xval),dim(xmad))
	expect_equal(length(xval),length(xmad))

	# do not error out
	show(xmad)

	set.char.seed("d35d4e4a-af3a-4491-a759-377fca599ec5")
	ddd <- matrix(rnorm(length(xval)*5),ncol=5)
	yt <- 'any'
	xt <- 'anx'
	vvv <- crossprod(matrix(rnorm(100*ncol(ddd)),ncol=ncol(ddd)))
	xmad <- madness(xval,ddd,vtag=yt,xtag=xt,varx=vvv)

	expect_equal(xval,val(xmad))
	expect_equal(ddd,dvdx(xmad))
	expect_equal(vvv,varx(xmad))
	expect_equal(xt,xtag(xmad))
	expect_equal(yt,vtag(xmad))

	# as.foo
	expect_equal(as.numeric(xval),as.numeric(xmad))
	expect_equal(as.matrix(xval),as.matrix(xmad))
	expect_equal(as.array(xval),as.array(xmad))

	expect_equal(as.array(xval),as(xmad,'array'))
	expect_equal(as.matrix(xval),as(xmad,'matrix'))
	expect_equal(as.numeric(xval),as(xmad,'numeric'))
	expect_equal(as.integer(xval),as(xmad,'integer'))
	expect_equal(as.logical(xval),as(xmad,'logical'))
	expect_equal(as.complex(xval),as(xmad,'complex'))

	#expect_equal(as.integer(xval),as.integer(xmad))
	#expect_equal(as.logical(xval),as.logical(xmad))
	#expect_equal(as.complex(xval),as.complex(xmad))

	ddd <- matrix(rnorm(length(xval)*5),ncol=5)
	yt <- 'anewy'
	xt <- 'anewx'
	vvv <- crossprod(matrix(rnorm(100*ncol(ddd)),ncol=ncol(ddd)))
	xtag(xmad) <- xt
	vtag(xmad) <- yt
	varx(xmad) <- vvv
	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("just vcov",{#FOLDUP
	yt <- 'any'
	xt <- 'anx'

	set.char.seed("b41fdb74-7fa4-4cbc-908d-4c7763403212")
	for (nr in c(1,4,8)) {
		xval <- matrix(1 + runif(nr*nr),nrow=nr)
		ddd <- matrix(rnorm(length(xval)*5),ncol=5)
		vvv <- crossprod(matrix(rnorm(100*ncol(ddd)),ncol=ncol(ddd)))
		xmad <- madness(xval,ddd,vtag=yt,xtag=xt,varx=vvv)
		blah <- vcov(xmad)
	}

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("just blockrep",{#FOLDUP
	yt <- 'any'
	xt <- 'anx'

	set.char.seed("b41fdb74-7fa4-4cbc-908d-4c7763403212")
	for (nr in c(1,4,8)) {
		xval <- matrix(1 + runif(nr*nr),nrow=nr)
		ddd <- matrix(rnorm(length(xval)*5),ncol=5)
		vvv <- crossprod(matrix(rnorm(100*ncol(ddd)),ncol=ncol(ddd)))
		xmad <- madness(xval,ddd,vtag=yt,xtag=xt,varx=vvv)
		blah <- blockrep(xmad,c(1))
		blah <- blockrep(xmad,c(1,1,1))
		blah <- blockrep(xmad,c(1,2,1))
		blah <- blockrep(xmad,c(2,1))
		blah <- blockrep(xmad,c(2,2,2))
		blah <- blockrep(xmad,c(2,2,1,2,3))
		expect_error(dumb <- blockrep(xmad,c(-1,1,1)))
		expect_error(dumb <- blockrep(xmad,c(1.5,1,1)))
		# xmad is nr x nr
		blah <- repto(xmad,c(nr,nr))
		blah <- repto(xmad,c(nr,nr,nr))
		blah <- repto(xmad,c(nr,2*nr))
		blah <- repto(xmad,c(nr,nr,3))
		blah <- repto(xmad,c(2*nr,3*nr,1,1,3))
		blah <- repto(xmad,c(4*nr,2*nr,1,3))
		expect_error(dumb <- repto(xmad,c(0.5*nr,-1,0)))
		expect_error(dumb <- repto(xmad,c(nr,0.5,1)))
		if (nr > 1) {
			expect_error(dumb <- repto(xmad,c(1,1,1)))
		}
	}

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("sumprodmaxmin?",{#FOLDUP
	yt <- 'any'
	xt <- 'anx'

	set.char.seed("820d69e3-9140-40b3-8b78-130d48a84b99")
	for (nr in c(4,6)) {
		xval <- matrix(1 + runif(nr*nr),nrow=nr)
		yval <- matrix(1 + runif(nr*nr),nrow=nr)
		xddd <- matrix(rnorm(length(xval)*5),ncol=5)
		yddd <- matrix(rnorm(length(yval)*5),ncol=5)

		for (pnan in c(0,25)) {
			if (pnan > 0) {
				xval[xval < 1 + (pnan/100)] <- NA
				yval[yval < 1 + (pnan/100)] <- NA
			}
			xmad <- madness(xval,xddd,vtag=yt,xtag=xt)
			ymad <- madness(yval,yddd,vtag=yt,xtag=xt)

			for (narm in c(FALSE,TRUE)) {
				blah <- sum(xmad,na.rm=narm)
				blah <- sum(xmad,ymad,na.rm=narm)
				blah <- sum(xmad,ymad,xval,na.rm=narm)
				blah <- prod(xmad,na.rm=narm)
				blah <- prod(xmad,ymad,na.rm=narm)
				blah <- prod(xmad,ymad,xval,na.rm=narm)

				blah <- max(xmad,na.rm=narm)
				blah <- max(xmad,ymad,na.rm=narm)
				blah <- max(ymad,xmad,na.rm=narm)
				blah <- max(xmad,ymad,xval,na.rm=narm)
				blah <- max(xmad,ymad,xval,yval,na.rm=narm)
				# these are *not* supported yet. bleah.
				#blah <- max(xmad,xval,ymad,na.rm=narm)
				#blah <- max(xval,xmad,ymad,na.rm=narm)

				blah <- min(xmad,na.rm=narm)
				blah <- min(xmad,ymad,na.rm=narm)
				blah <- min(ymad,xmad,na.rm=narm)
				blah <- min(xmad,ymad,xval,na.rm=narm)
				blah <- min(xmad,ymad,xval,yval,na.rm=narm)
				# these are *not* supported yet. bleah.
				#blah <- min(xmad,xval,ymad,na.rm=narm)
				#blah <- min(xval,xmad,ymad,na.rm=narm)
			}
		}
	}

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("outer?",{#FOLDUP
	yt <- 'any'
	xt <- 'anx'

	set.char.seed("b41fdb74-7fa4-4cbc-908d-4c7763403212")
	for (nr in c(1,4,8)) {
		xval <- matrix(1 + runif(nr*nr),nrow=nr)
		ddd <- matrix(rnorm(length(xval)*5),ncol=5)
		vvv <- crossprod(matrix(rnorm(100*ncol(ddd)),ncol=ncol(ddd)))
		xmad <- madness(xval,ddd,vtag=yt,xtag=xt,varx=vvv)

		yval <- matrix(1 + runif(nr*3),nrow=nr)
		ddd <- matrix(rnorm(length(yval)*5),ncol=5)
		vvv <- crossprod(matrix(rnorm(100*ncol(ddd)),ncol=ncol(ddd)))
		ymad <- madness(yval,ddd,vtag=yt,xtag=xt,varx=vvv)

		blah <- outer(xmad,ymad,'*')
		blah <- outer(xmad,ymad,'+')
		blah <- outer(xmad,ymad,'-')
	}

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("vech",{#FOLDUP
	# first on arrays...
	xv <- array(rnorm(3),dim=c(3,1))
	MV <- array(rnorm(9),dim=c(3,3))


	vec(MV)
	vech(MV)
	vech(MV,1)
	vech(MV,-1)

	ivech(xv)
	ivech(xv,symmetric=TRUE)
	ivech(xv,-1)
	ivech(xv,-1,symmetric=TRUE)

	ev <- array(rnorm(8),dim=c(8,1))
	ivech(ev,1)
	vech(ivech(ev,1),1)

	set.char.seed("b41fdb74-7fa4-4cbc-908d-4c7763403212")
	for (nr in c(3,6,10)) {
		yt <- 'any'
		xt <- 'anx'
		xval <- matrix(1 + runif(nr),nrow=nr)
		ddd <- matrix(rnorm(length(xval)*5),ncol=5)
		vvv <- crossprod(matrix(rnorm(100*ncol(ddd)),ncol=ncol(ddd)))
		xmad <- madness(xval,ddd,vtag=yt,xtag=xt,varx=vvv)

		blah <- ivech(xmad)
		blah <- ivech(xmad,-1)
		blah <- ivech(xmad,-1,symmetric=TRUE)
	}

	for (nr in c(8,13)) {
		yt <- 'any'
		xt <- 'anx'
		xval <- matrix(1 + runif(nr),nrow=nr)
		ddd <- matrix(rnorm(length(xval)*5),ncol=5)
		vvv <- crossprod(matrix(rnorm(100*ncol(ddd)),ncol=ncol(ddd)))
		xmad <- madness(xval,ddd,vtag=yt,xtag=xt,varx=vvv)

		blah <- ivech(xmad,1)
		blah <- ivech(xmad,1,symmetric=FALSE)
		# cannot symmetric when k > 0
		expect_error(blah <- ivech(xmad,1,symmetric=TRUE))
	}

	# these should error out! wrong size!
	for (nr in c(16,27)) {
		yt <- 'any'
		xt <- 'anx'
		xval <- matrix(1 + runif(nr),nrow=nr)
		ddd <- matrix(rnorm(length(xval)*5),ncol=5)
		vvv <- crossprod(matrix(rnorm(100*ncol(ddd)),ncol=ncol(ddd)))
		xmad <- madness(xval,ddd,vtag=yt,xtag=xt,varx=vvv)

		expect_error(dumb <- ivech(xmad,0))
		expect_error(dumb <- ivech(xmad,-1))
	}

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("eigen?",{#FOLDUP
	yt <- 'any'
	xt <- 'anx'
	set.char.seed("b41fdb74-7fa4-4cbc-908d-4c7763403212")

	for (nr in c(4)) {
		xval <- matrix(1 + runif(nr*nr),nrow=nr)
		ddd <- matrix(rnorm(length(xval)*5),ncol=5)
		vvv <- crossprod(matrix(rnorm(100*ncol(ddd)),ncol=ncol(ddd)))
		xmad <- madness(xval,ddd,vtag=yt,xtag=xt,varx=vvv)
		foo <- eigen(xmad + t(xmad))
	}

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("theta",{#FOLDUP
	# first on arrays...
	set.char.seed('d8ccbb36-1002-4c9e-81d9-0ee6b173047a')
	MV <- array(rnorm(100*3),dim=c(100,3))
	th <- theta(MV)
	MV <- array(rnorm(100*3*3),dim=c(100,3,3))
	th <- theta(MV)
	MV <- array(rnorm(100*3*3*3),dim=c(100,3,3,3))
	th <- theta(MV)
	th <- theta(MV,xtag='FOO')

	# and as data frame.
	MV <- data.frame(a=runif(100),b=rnorm(100),c=exp(runif(100)))
	th <- theta(MV)
	th <- theta(MV,xtag='FOO')

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("twomoments",{#FOLDUP
	# first on arrays...
	set.char.seed('34b7717b-ddff-473a-98d1-422b040b6d82')

	checkem <- list(array(rnorm(100*3),dim=c(100,3)),
								 	array(rnorm(100*3*3),dim=c(100,3,3)),
									array(rnorm(100*3*3*3),dim=c(100,3,3,3)))

	for (ccc in seq_along(checkem)) {
		MV <- checkem[[ccc]]
		th <- twomoments(MV)
		th <- twomoments(MV,xtag='FOO')
		th <- twomoments(MV,xtag='FOO',df=0)
		th <- twomoments(MV,xtag='FOO',df=0,diag.only=TRUE)
		th <- twomoments(MV,xtag='FOO',df=0,diag.only=TRUE,vcov.func=vcov)
	}

	# and as data frame.
	MV <- data.frame(a=runif(100),b=rnorm(100),c=exp(runif(100)))
	th <- twomoments(MV)
	th <- twomoments(MV,xtag='FOO')
	th <- twomoments(MV,xtag='FOO',df=0)
	th <- twomoments(MV,xtag='FOO',df=0,diag.only=TRUE)

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("to_objective",{#FOLDUP
	set.char.seed('1407dc0a-ca1b-4620-b436-e13df111b6ab')
	
	MV <- array(rnorm(100*3),dim=c(100,3))
	madM <- madness(MV)
	Mnorm <- norm(madM)
	to_objective(Mnorm)
	expect_error(dumb <- to_objective(madM))

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
test_that("numerical derivative",{#FOLDUP
	# first on arrays...
	set.char.seed('818fae1e-30b3-4bcb-bb54-691f4c8d05ae')
	MV <- array(rnorm(100*3),dim=c(100,3))
	madv <- madness(MV)
	fn <- function(x) { cos(x) }
	gn <- function(x,h) { exp(x + h) }
	for (ty in c('forward','backward','central')) {
		for (eps in c(1e-4,1e-9)) {
			dMV1 <- numderiv(fn,MV,type=ty,eps=eps)
			dMV2 <- numderiv(gn,MV,type=ty,eps=eps,h=-3)
			dMV3 <- numderiv(fn,madv,type=ty,eps=eps)
			dMV4 <- numderiv(gn,madv,type=ty,eps=eps,h=-3)
		}
	}

	# sentinel:
	expect_true(TRUE)
})#UNFOLD

#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
