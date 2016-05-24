# Unit tets for RNGseq
# 
# Author: Renaud Gaujoux
###############################################################################

library(parallel)

test.RNGseq_seed <- function(){
	
	# actual testing function
	.test_loc <- function(.msg, ..., .change=FALSE){
		msg <- function(...) paste(.msg, ':', ...)
		os <- RNGseed()
		on.exit(RNGseed(os))
		s <- RNGseq_seed(...)
		checkTrue(length(s) == 7L && s[1] %% 100 == 7L, msg("RNGseq_seed returns a value of .Random.seed for L'Ecuyer-CMRG"))
		checkIdentical(RNGseed()[1], os[1], msg("RNGseq_seed does not change the type of RNG"))
		
		if( !.change ) checkIdentical(RNGseed(), os, msg("RNGseq_seed does not change the value of .Random.seed"))
		else checkTrue( !identical(RNGseed(), os), msg("RNGseq_seed changes the value of .Random.seed"))
		s
	}
	
	# test in two RNG settings: default and L'Ecuyer
	.test <- function(.msg, ..., ss=NULL, .change=FALSE, Dchange=.change, Lchange=.change){
		os <- RNGseed()
		on.exit(RNGseed(os))
		
		# default RNG
		RNGkind('default')
		if( !is.null(ss) ) set.seed(ss)
		s1 <- .test_loc(paste(.msg, '- default'), ..., .change=Dchange)
		
		RNGkind("L'Ecuyer")
		if( !is.null(ss) ) set.seed(ss)
		s2 <- .test_loc(paste(.msg, "- CMRG"), ..., .change=Lchange)
		
		list(s1, s2)
	}
	
	os <- RNGseed()
	on.exit(RNGseed(os))
	
	RNGkind('default', 'default')
	
	# test different arguments
	s1 <- .test("seed=missing", ss=1, Dchange=TRUE, Lchange=FALSE)
	runif(10)
	s2 <- .test("seed=NULL", NULL, ss=1, Dchange=TRUE, Lchange=FALSE)
	checkIdentical(s1, s2, "set.seed(1) + seed=missing and seed=NULL return identical results")
	
	# doRNG seed with single numeric
	runif(10)
	s3 <- .test("seed=single numeric", 1)
	checkIdentical(s1[[1]], s3[[1]], "v1.4 - set.seed(1) + seed=missing and seed=1 return identical results when current RNG is NOT CMRG")
	checkIdentical(s1[[2]], s3[[2]], "v1.4 - set.seed(1) + seed=missing and seed=1 return identical results when current RNG is CMRG")
	checkTrue( !identical(s1[[1]], s1[[2]]), "v1.4 - set.seed(1) + seed=missing return NON identical results in different RNG settings")
	checkTrue( !identical(s3[[1]], s3[[2]]), "v1.4 - seed=num return NON identical results in different RNG settings")
	
	# version < 1.4
#	doRNGversion("1.3.9999")
	s1 <- .test("v1.3 - seed=missing", ss=1, Dchange=TRUE, Lchange=TRUE, version=1)
	s3 <- .test("v1.3 - seed=single numeric", 1, version=1)
	checkIdentical(s1[[1]], s3[[1]], "v1.3 - set.seed(1) + seed=missing and seed=1 return identical results when current RNG is NOT CMRG")
	checkTrue( !identical(s1[[2]], s3[[2]]), "v1.3 - set.seed(1) + seed=missing and seed=1 return NON identical results when current RNG is CMRG")
	checkTrue( !identical(s1[[1]], s1[[2]]), "v1.3 - set.seed(1) + seed=missing return NON identical results in different RNG settings")
	checkTrue( !identical(s3[[1]], s3[[2]]), "v1.4 - seed=num return NON identical results in different RNG settings")
#	doRNGversion(NULL) 
	##
	
	.test("seed=single integer", 10L)
	# directly set doRNG seed with a 6-length
	.test("seed=6-length integer", 1:6)
	.test("seed=6-length numeric", as.numeric(1:6))
	s <- 1:6
	checkIdentical(RNGseq_seed(s)[2:7], s, "RNGseq_seed(6-length) returns stream to the given value")
	# directly set doRNG seed with a full 7-length .Random.seed
	.test("seed=7-length integer", c(407L,1:6))
	.test("seed=7-length numeric", as.numeric(c(107L,1:6)))
	s <- c(407L,1:6)
	checkIdentical(RNGseq_seed(s), s, "RNGseq_seed(7-length) returns complete seed with the given value")
	
	# errors
	os <- RNGseed()
	checkException(RNGseq_seed(NA), "seed=NA throws an exception")
	checkIdentical(os, RNGseed(), "RNGseq_seed(NA) does not change the value of .Random.seed [error]")
	
	# Current CMRG is L'Ecuyer
	RNGkind("L'Ecuyer")
	set.seed(456)
	s <- RNGseed()
	r <- RNGseq_seed(NULL)
	checkIdentical(s, r, "Current is CMRG: seed=NULL return current stream")
	runif(10)
	checkIdentical(s, RNGseq_seed(456), "Current is CMRG: seed=numeric return stream seeded with value")
	
}

test.RNGseq <- function(){
	
	os <- RNGseed()
	on.exit(RNGseed(os))
	
	# actual testing function
	.test_loc <- function(.msg, n, ..., .list=TRUE, .change=FALSE){
		msg <- function(...) paste(.msg, ':', ...)
		os <- RNGseed()
		on.exit(RNGseed(os))
	
		s <- RNGseq(n, ...)
		
		if( !.change ) checkIdentical(RNGseed(), os, msg("the value of .Random.seed is not changed"))
		else checkTrue( !identical(RNGseed(), os), msg("the value of .Random.seed does change"))
		
		if( .list )	checkTrue(is.list(s), msg("result is a list"))
		else{
			checkTrue(is.integer(s), msg("result is an integer vector"))
			s <- list(s)
		}
		
		checkTrue(length(s) == n, msg("result has correct length"))
		checkTrue(all(sapply(s, length) == 7L), msg("each element has length 7"))
		checkTrue(all(sapply(s, function(x) x[1] %% 100) == 7L), msg("each element has correct RNG kind"))
		s
	}
	
	.test <- function(msg, n, ...){
		set.seed(1)
		s1 <- .test_loc(paste(msg, '- no seed'), n, ..., .change=TRUE)
		runif(1)
		s2 <- .test_loc(paste(msg, '- seed=1'), n, 1, ..., .change=FALSE)
		#checkIdentical(s1, s2, paste(msg, " - set.seed(1) + no seed is identical to seed=1"))
		.test_loc(paste(msg, '- seed=1:6'), n, 1:6, ...)
	}
	.test("n=1", 1, .list=FALSE)
	.test("n=2", 2)
	.test("n=5", 5)
	
	# with full list
	s <- RNGseq(3)
	checkIdentical(RNGseq(length(s), s), s, "If passing a complete list: returns the list itself")
	s3 <- RNGseq(5)
	s <- structure(s, rng=s3)
	checkIdentical(RNGseq(length(s3), s), s3, "If passing a complete list in rng S3 slot: returns the complete slot")
	#

	# Current RNG is CMRG
	set.seed(456, "L'Ec")
	s <- .Random.seed
	ref <- list(s, nextRNGStream(s), nextRNGStream(nextRNGStream(s)))
	rs <- RNGseq(3, 456)
	checkIdentical(rs, ref, "Current RNG is CMRG: RNGseq(n, num) returns RNG streams that start with stream as set.seed")
	checkIdentical(s, .Random.seed, "Current RNG is CMRG: RNGseq(n, num) did not change random seed")
	
	runif(10)
	s <- .Random.seed
	ref <- list(s, nextRNGStream(s), nextRNGStream(nextRNGStream(s)))
	rs2 <- RNGseq(3)
	checkIdentical(rs2, ref, "Current RNG is CMRG: RNGseq(n) returns RNG streams that start with current stream")
	checkIdentical(.Random.seed, nextRNGStream(tail(rs2,1)[[1]]), "Current RNG is CMRG: RNGseq(n) changes current random seed to next stream of last stream in sequence")
	
}
