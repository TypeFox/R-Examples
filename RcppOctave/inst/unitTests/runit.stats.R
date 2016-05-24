# Unit tests for stat functions.
#
# Copyright (C) 2011 Renaud Gaujoux
# 
# This file is part of RcppOctave.
#
# RcppOctave is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RcppOctave is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RcppOctave.  If not, see <http://www.gnu.org/licenses/>.


# Auxiliary function for testing the o_r* functions
check.rfun <- function(rfun, ofun, ..., .msg=NULL){
	
	rname <- deparse(substitute(rfun))
	oname <- deparse(substitute(ofun))
	title <- paste(rname, ' and ', oname, if( !missing(.msg) ) paste('-', .msg),  ':', sep='')
	msg <- function(...) paste(title, ...)
	
	checkRand <- function(e1, e2, msg){
		set.seed(123)
		irs1 <- .Random.seed
		r1 <- force(e1)
		rs1 <- .Random.seed
		
		set.seed(123);
		irs2 <- .Random.seed
		r2 <- force(e2)
		rs2 <- .Random.seed
		
		checkIdentical(irs1, irs2, paste(msg, ": .Random.seed identical after set.seed, before draws"))
		checkTrue( !identical(irs2, rs2), paste(msg, ": .Random.seed changed after Octave draw"))
		checkIdentical(r1, r2, paste(msg, ": result is identical"))
		checkIdentical(rs1, rs2, paste(msg, ": .Random.seed identical after Octave draw and R draw"))
	}
	
	checkRand( rfun(10, ...), ofun(1,10, ...), msg('random vector'))
	checkRand( matrix(rfun(25, ...), 5, 5), ofun(5, ...), msg('random square matrix'))
	checkRand( matrix(rfun(30, ...), 10, 3), ofun(10, 3, ...), msg('random rectangular matrix'))
}

#' Common unit test for o_runif
test.o_runif <- function(){
	check.rfun(runif, .O$rand)
	check.rfun(runif, o_runif)
}

#' Unit test for o_rnorm
test.o_rnorm <- function(){	
	check.rfun(rnorm, .O$randn)
	check.rfun(rnorm, o_rnorm)
}

#' Unit test for o_rexp
test.o_rexp <- function(){	
	check.rfun(rexp, .O$rande)
	check.rfun(rexp, o_rexp)
}

#' Unit test for o_rgamma
test.o_rgamma <- function(){
	check.rfun(rgamma, o_rgamma, shape=1, scale=1, .msg='shape=1 / scale=1')
	check.rfun(rgamma, o_rgamma, shape=10, scale=10, .msg='shape=10 / scale=10')
	check.rfun(rgamma, o_rgamma, shape=10, scale=1, .msg='shape=10 / scale=1')
	check.rfun(rgamma, o_rgamma, shape=1, scale=10, .msg='shape=1 / scale=10')
}

#' Unit test for o_rpois
test.o_rpois <- function(){
	check.rfun(rpois, o_rpois, lambda=4, .msg='lambda=4')
}
