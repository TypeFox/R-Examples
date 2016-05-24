# Unit test for doRNG
# 
# Author: Renaud Gaujoux
# Creation: 28 Mar 2012
###############################################################################


#test.CMRGseed <- function(){
#	
#	msg <- function(...) paste(.msg, ':', ...)
#	
#	# Unit tests
#	.msg <- "Call CMRGseed without argument"
#	rs <- .Random.seed
#	checkIdentical( length(CMRGseed()), 7L, msg("Seed is of length 7") )
#	checkIdentical(rs, .Random.seed, msg("does not change .Random.seed"))
#	
#	.msg <- "Call CMRGseed with a single argument"
#	rs <- .Random.seed
#	checkIdentical( length(CMRGseed(1)), 7L, msg("Seed is of length 7") )
#	checkIdentical(rs, .Random.seed, msg("does not change .Random.seed"))
#	checkTrue( all(!is.na(CMRGseed(1))), msg("No NA in the returned seed") )
#	
#}



test.dorng <- function(){
	
		test_dopar <- function(.msg, s.seq){
			
			orng <- RNGseed()
			on.exit({ doRNGversion(NULL); RNGseed(orng); registerDoSEQ()} )
			
			msg <- function(...) paste(.msg, ':', ...)
			noattr <- function(x){ attributes(x) <- NULL; x}
			
			# standard %dopar% loops are _not_ reproducible
			set.seed(1234)
			s1 <- foreach(i=1:4) %dopar% { runif(1) }
			set.seed(1234)
			s2 <- foreach(i=1:4) %dopar% { runif(1) }
			if( !missing(s.seq) ) checkTrue( !identical(s1, s2), msg("Standard %dopar% loop is not reproducible"))
			
			# %dorng% loops ensure reproducibility
			local({
				set.seed(1234)
				s1 <- foreach(i=1:4) %dorng% { runif(1) }
				runif(10)
				set.seed(1234)
				s2 <- foreach(i=1:4) %dorng% { runif(1) }
				checkIdentical(s1, s2, msg("%dorng% loop is reproducible with set.seed"))
				# check RNG settings in result
				set.seed(1234)
				ref <- RNGseq(4)
				rngs <- attr(s1, 'rng')
				checkTrue(!is.null(rngs), msg("Results contains RNG data"))
				checkIdentical(rngs, ref, msg("Results contains whole sequence of RNG seeds"))
			})
			
			# or 
			local({
				s1 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
				s2 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
				checkIdentical(s1, s2, msg("%dorng% loop is reproducible with .options.RNG"))
				ref <- RNGseq(4, 1234)
				rngs <- attr(s1, 'rng')
				checkTrue(!is.null(rngs), msg("Results of seed loop contains RNG data"))
				checkIdentical(rngs, ref, msg("Results of seeded loop contains whole sequence of RNG seeds"))
			})
			
			## check extra arguments to .options.RNG
			# Normal RNG parameter is taken into account
			s.unif.noNk <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(5) }
			s.unif.wNk <- foreach(i=1:4, .options.RNG=list(1234, normal.kind="Ahrens-Dieter")) %dorng% { runif(5) }
			checkIdentical(noattr(s.unif.noNk), noattr(s.unif.wNk)
				, msg("%dorng% loop (runif) with normal.kind in .options.RNG is identical as without"))
		
			s.norm.noNk <- foreach(i=1:4, .options.RNG=1234) %dorng% { rnorm(5) }
			s.norm.wNk <- foreach(i=1:4, .options.RNG=list(1234, normal.kind="Ahrens-Dieter")) %dorng% { rnorm(5) }
			checkTrue(!isTRUE(all.equal(noattr(s.norm.noNk), noattr(s.norm.wNk)))
				, msg("%dorng% loop (rnorm) with normal.kind in .options.RNG is different as without"))			

			# reproduce previous loop
			runif(10)
			s1 <- foreach(i=1:4) %dorng% { runif(5) }
			s2 <- foreach(i=1:4, .options.RNG=s1) %dorng% { runif(5) }
			checkIdentical(s1, s2, "Seeding using .options.RNG={result from other loop} give identical results")
			# directly set sequence of seeds
			sL <- list(c(407,1:6), c(407,11:16), c(407,21:26), c(407, 31:36))
			s1.L <- foreach(i=1:4, .options.RNG=sL) %dorng% { runif(5) }
			runif(10)
			s2.L <- foreach(i=1:4, .options.RNG=sL) %dorng% { runif(5) }
			checkIdentical(s1.L, s2.L, "Seeding using .options.RNG=list twice give identical results")
			# directly set sequence of seeds as a matrix
			sM <- sapply(sL, identity)
			s1.M <- foreach(i=1:4, .options.RNG=sM) %dorng% { runif(5) }
			runif(10)
			s2.M <- foreach(i=1:4, .options.RNG=sM) %dorng% { runif(5) }
			checkIdentical(s1.M, s2.M, "Seeding using .options.RNG=matrix twice give identical results")
			checkIdentical(s1.M, s1.L, "Seeding using .options.RNG=matrix gives identical results as the same seed in list")
			
			# separate %dorng% loops are different
			set.seed(1234)
			rs1 <- .Random.seed
			s1 <- foreach(i=1:4) %dorng% { runif(1) }
			rs1_2 <- .Random.seed
			s2 <- foreach(i=1:4) %dorng% { runif(1) }
			checkTrue( !identical(rs1, rs1_2), msg("unseed %dorng% loop changes .Random.seed"))
			checkTrue( !identical(s1, s2), msg("two consecutive (unseeded) %dorng% loops are not identical"))
			checkTrue( !identical(unlist(s1), unlist(s2)), msg("two consecutive (unseeded) %dorng% loops are not identical (values)"))
			
			# But the whole sequence of loops is reproducible
			set.seed(1234)
			s1.2 <- foreach(i=1:4) %dorng% { runif(1) }
			s2.2 <- foreach(i=1:4) %dorng% { runif(1) }
			checkTrue( identical(s1, s1.2) && identical(s2, s2.2), msg("set.seed + two consecutive %dorng% loops are reproducible"))
			
			s <- list(s1, s2)
			if( !missing(s.seq) )
				checkIdentical(s, s.seq, msg("result is identical to sequential computation"))
			
			# check behaviour with set.seed
			set.seed(789)
			s1 <- foreach(i=1:6) %dorng%{ runif(1) }
			s2 <- foreach(i=1:6, .options.RNG=789) %dorng%{ runif(1) }
			checkIdentical(s1, s2, "set.seed before %dorng% is identical to using .options.RNG")
			# current RNG is CRMG
			set.seed(789, "L'Ec")
			s1 <- foreach(i=1:6) %dorng%{ runif(1) }
			s2 <- foreach(i=1:6, .options.RNG=789) %dorng%{ runif(1) }
			checkIdentical(s1, s2, "set.seed before %dorng% is identical to using .options.RNG, if current RNG is CRMG")

			s
		}
		
		orng <- RNGseed()
		on.exit({ doRNGversion(NULL); RNGseed(orng); registerDoSEQ()} )
		
		library(doParallel)
		
		# Sequential computation
		registerDoSEQ()
		s.seq <- test_dopar("Sequential")
		
		# Multicore cluster
    if( .Platform$OS.type != 'windows'){
		  registerDoParallel(cores=2)
		  s <- test_dopar("Multicore", s.seq)
    }
		
		# SNOW-like cluster
		cl <- makeCluster(2)
		registerDoParallel(cl)
		test_dopar("SNOW-like cluster", s.seq)
		stopCluster(cl)
		
		# Works with doMPI
		if( require(doMPI) ){
			cl <- startMPIcluster(2)
			registerDoMPI(cl)
			test_dopar("MPI cluster", s.seq) 
			closeCluster(cl)
		}
}

test.registerDoRNG <- function(){
	
	on.exit( stopCluster(cl) )
	library(doParallel)
	cl <- makeCluster(2)
	registerDoParallel(cl)
	
	# One can make existing %dopar% loops reproducible using %dorng% loops or registerDoRNG  
	set.seed(1234)
	r1 <- foreach(i=1:4) %dorng% { runif(1) }
	registerDoRNG()
	set.seed(1234)
	r2 <- foreach(i=1:4) %dopar% { runif(1) }
	checkIdentical(r1, r2, "registerDoRNG + set.seed: makes a %dopar% loop behave like a set.seed + %dorng% loop")
	stopCluster(cl)
	
	# Registering another foreach backend disables doRNG
	cl <- makeCluster(2)
	registerDoParallel(cl)
	set.seed(1234)
	s1 <- foreach(i=1:4) %dopar% { runif(1) }
	set.seed(1234)
	s2 <- foreach(i=1:4) %dorng% { runif(1) }
	checkTrue( !identical(s1, s2), "Registering another foreach backend disables doRNG")
	
	# doRNG is re-nabled by re-registering it 
	registerDoRNG()
	set.seed(1234)
	r3 <- foreach(i=1:4) %dopar% { runif(1) }
	checkIdentical(r2, r3, "doRNG is re-nabled by re-registering it")
	r4 <- foreach(i=1:4) %dopar% { runif(1) }
	# NB: the results are identical independently of the task scheduling
	# (r2 used 2 nodes, while r3 used 3 nodes)
	
	# Reproducibility of sequences of loops
	# pass seed to registerDoRNG 
	runif(10)
	registerDoRNG(1234)
	r1 <- foreach(i=1:4) %dopar% { runif(1) }
	r2 <- foreach(i=1:4) %dopar% { runif(1) }
	registerDoRNG(1234)
	r3 <- foreach(i=1:4) %dopar% { runif(1) }
	r4 <- foreach(i=1:4) %dopar% { runif(1) }
	checkIdentical(r3, r1, "registerDoRNG(1234) allow reproducing sequences of %dopar% loops (1)")
	checkIdentical(r4, r2, "registerDoRNG(1234) allow reproducing sequences of %dopar% loops (2)")	
	# use set.seed
	runif(10)
	registerDoRNG()
	set.seed(1234)
	s1 <- foreach(i=1:4) %dopar% { runif(1) }
	s2 <- foreach(i=1:4) %dopar% { runif(1) }
	set.seed(1234)
	s3 <- foreach(i=1:4) %dopar% { runif(1) }
	s4 <- foreach(i=1:4) %dopar% { runif(1) }
	checkIdentical(s3, s1, "registerDoRNG + set.seed(1234) allow reproducing sequences of %dopar% loops (1)")
	checkIdentical(s4, s2, "registerDoRNG + set.seed(1234) allow reproducing sequences of %dopar% loops (2)")
	
	runif(5)
	registerDoRNG()
	set.seed(1234)
	s5 <- foreach(i=1:4) %dopar% { runif(1) }
	s6 <- foreach(i=1:4) %dopar% { runif(1) }
	checkIdentical(s5, r3, "registerDoRNG() + set.seed give same results as registerDoRNG(1234) (1)")
	checkIdentical(s6, r4, "registerDoRNG() + set.seed give same results as registerDoRNG(1234) (2)")
	
	# argument `once=FALSE` reseed doRNG's seed at the beginning of each loop 
	registerDoRNG(1234, once=FALSE)
	r1 <- foreach(i=1:4) %dopar% { runif(1) }
	r2 <- foreach(i=1:4) %dopar% { runif(1) }
	r3 <- foreach(i=1:4, .options.RNG=1234) %dopar% { runif(1) }
	checkIdentical(r1, r2, "argument `once=FALSE` reseed doRNG's seed at the beginning of each loop")
	checkIdentical(r1, r3, "argument `once=FALSE` reseed %dorng% loop as .options.RNG")
	
	# Once doRNG is registered the seed can also be passed as an option to %dopar%
	r1.2 <- foreach(i=1:4, .options.RNG=456) %dopar% { runif(1) }
	r2.2 <- foreach(i=1:4, .options.RNG=456) %dopar% { runif(1) }	
	checkIdentical(r1.2, r2.2, "Once doRNG is registered the seed can also be passed as an option to %dopar%")
	checkTrue(!identical(r1.2, r1), "The seed passed as an option is really taken into account")
	
}
