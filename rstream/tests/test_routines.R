#############################################################################
##                                                                         ##
##   Test routines used for all rstream classes                            ##
##                                                                         ##
#############################################################################

## Remark -------------------------------------------------------------------

## This file only defines test routines but does not execute any test.
## It should be loaded by the test files:
##    source ("test_routines.R")
## The default test parameters can be overwritten.


## Test Parameters ----------------------------------------------------------

## API
if (!exists("samplesize")) samplesize <- 1e4
if (!exists("streamname")) streamname <- "testname"

## Goodness-of-Fit tests
if (!exists("chi2.samplesize"))  chi2.samplesize <- 1e6
if (!exists("alpha"))            alpha <- 1.e-3


## Load library -------------------------------------------------------------

library(rstream)


## Test API of package ------------------------------------------------------

rstream.check.API <- function (
                               type,		  # type (class) of stream 
                               ...,		  # optional args for new(...)
                               t_sample=TRUE,
                               t_reset=TRUE,
                               t_clone=TRUE,
                               t_pack=TRUE,
                               t_name=TRUE, 
                               t_anti=TRUE,
                               t_incp=TRUE, 
                               t_print=TRUE ) {

	## Create a stream ..................................................
	s <- new (type, name=streamname, ...)


	## Print ............................................................
	print(s)
        
        
	## Sample ...........................................................
	x <- rstream.sample(s,samplesize)
	if (length(x) != samplesize)
		stop("sample.rstream failed:1")
	
        
	## Reset ............................................................
	rstream.reset(s)
	y <- rstream.sample(s,samplesize)
	if (!identical(all.equal(x, y), TRUE))
		stop("reset failed:1")
        
        
	## Sample: r ........................................................
	rstream.reset(s)
	y <- r(s,samplesize)
	if (!identical(all.equal(x, y), TRUE))
		stop("r failed:1")
        
        
	## Name ............................................................
	if (streamname != rstream.name(s))
		stop("name failed:1")
        
	xyzname <- paste("xyz","streamname",sep="")
	rstream.name(s) <- xyzname
	if (xyzname != rstream.name(s))
		stop("name failed:2")
        
	rstream.name(s) <- streamname  ## restore stream
        
        
	## Antithetic variables .............................................
	rstream.reset(s)
	state <- store_state(s)
        
	if(rstream.antithetic(s))
		stop("antithetic must be FALSE at default:1")
        
	rstream.reset(s)
	rstream.antithetic(s) <- TRUE
	if(!rstream.antithetic(s))
		stop("antithetic must be TRUE now:2")
	if(!comp_states(state, store_state(s), anti=TRUE))
		stop("antithetic failed:3")
        
	rstream.reset(s)
	rstream.antithetic(s) <- FALSE
	if(rstream.antithetic(s))
		stop("antithetic must be FALSE now:4")
	if(!comp_states(state, store_state(s)))
		stop("antithetic failed:5")
        
        
	## from now on rstream.antithetic(s) == TRUE
	rstream.antithetic(s) <- TRUE
	rstream.reset(s)
	state <- store_state(s)
        
        
	## Incprecision .....................................................
	if(rstream.incprecision(s))
		stop("incprecision must be FALSE at default:1")
	if (t_incp) {
		rstream.incprecision(s) <- TRUE
		if(!rstream.incprecision(s))
			stop("incprecision must be TRUE now:2")
		rstream.incprecision(s) <- FALSE
		if(rstream.incprecision(s))
                        stop("incprecision must be FALSE now:3") 
                
		rstream.reset(s)
		if(!comp_states(state, store_state(s)))
			stop("incprecision failed:4") }
        
        
	## Clone ............................................................
	rstream.reset(s)
	clone <- rstream.clone(s)
	if (rstream.name(clone) != paste(streamname,".",sep=""))
		stop("name of clone failed:1")
        
	if(!comp_states(state, store_state(clone),clone=TRUE))
                stop("cloning failed:2")
	
	rstream.reset(s)
	rstream.reset(clone)
	rstream.sample(clone,1)
	rstream.antithetic(clone) <- FALSE
	if (t_incp) rstream.incprecision(clone) <- TRUE
	if(!comp_states(state, store_state(s)))
		stop("cloned object not independent:3")
        
        
	## Pack/unpack ......................................................
	rstream.reset(s)
	rstream.packed(s)
	if(rstream.packed(s))
		stop("packed must be FALSE at default:1")
	rstream.packed(s) <- TRUE
	print(s)
	if(!rstream.packed(s))
		stop("packed must be TRUE now:2")
	dup <- s
	rstream.packed(s) <- FALSE
	if(!comp_states(state, store_state(s)))
		stop("packing failed:3")
        
	rstream.packed(dup) <- FALSE
	if(!comp_states(state, store_state(dup)))
		stop("packing failed:4")
	
	## "<-" .............................................................
	copy <- s
	rstream.reset(copy)
	if(!comp_states(state, store_state(copy)))
		stop("\"<-\" failed:1")

	rstream.reset(s)
	rstream.name(copy) <- "karl"
	rstream.antithetic(copy) <- FALSE
	if (t_incp) rstream.incprecision(copy) <- TRUE
	if(comp_states(state, store_state(s)))
		stop("\"<-\" failed:2")
        
	rstream.reset(copy)
	state1 <- store_state(s)
	rstream.reset(copy)
	if(!comp_states(state1, store_state(s)))
		stop("\"<-\" failed:3")
}


## Test interface rstream <-> R RNG -----------------------------------------

rstream.check.setRNG <- function (
                                  type,		  # type (class) of stream 
                                  ...) {	  # optional args for new(...)

        ## Create a stream ..................................................
	s <- new (type, name=streamname, ...)

        ## make a working copy
        sc <- rstream.clone(s)
        
	## set RNG
	rstream.RNG(s)
        
	## s und runif() should produce the same sequence
	x <- rstream.sample(sc, samplesize)
	y <- runif(samplesize)
	if (!identical(all.equal(x, y), TRUE))
		stop("set RNG failed: 1")

        ## try again
        ## (required as we had a bug in all pre-1.3 versions!)
	rstream.RNG(s)
	y <- runif(samplesize)
	if (!identical(all.equal(x, y), TRUE))
		stop("set RNG failed: 2")
        
	## get RNG
	sr <- rstream.RNG()
	x <- rstream.sample(sr, samplesize)
	rstream.antithetic(sr) <- TRUE   # this must not influence runif()
	y <- runif(samplesize)
	if (!identical(all.equal(x, y), TRUE))
		stop("set RNG failed: 3")
        
}


## Store state of generator in list -----------------------------------------

store_state <- function(stream) {
	## make empty list
	state <- list()
	## store data about stream
	state$name <- rstream.name(stream)
	state$anti <- rstream.antithetic(stream)
	state$incp <- rstream.incprecision(stream)
	## for getting information about the internal state we use a sample
	## (we cannot store these data directly since the structure of
	## generators is too different).
	## Moreover, we also want to check whether the generator works
	## correctly
	state$samp <- rstream.sample(stream, samplesize)
	## return state
	state
} # end of store_state()


## Compare state of generator (stored in list) ------------------------------

comp_states <- function(stateA, stateB, anti=FALSE, clone=FALSE) {
	## if 'anti' is TRUE then the antithetic flags in stateA and stateB
	## must be different. Otherwise they must coincide.
	if (anti) {
		stateB$anti <- !stateB$anti
		stateB$samp <- (1 - stateB$samp) }
        
        ## if 'clone' is TRUE then we must add a period '.' to the name.
	if (clone)
		stateA$name <- paste(stateA$name,".",sep="")
        
	## return result (TRUE if equal, FALSE if they differ)
	return (identical(all.equal(stateA, stateB), TRUE))
} # end of comp_states()


## --- Function for running chi^2 goodness-of-fit test ----------------------

rstream.check.chi2 <- function (type, ...) {
        ##  Run a chi^2 test and evaluate p-value.
        ##  
        ##  type   ...  type (class) of stream 
        ##  ...    ...  optional args for new(...)
        ##

        for (i in 1:2) {
  
          ## -- run test
          pval <- rstream.check.chi2.run(type, ...)

          ## -- check p-value
          if (pval > alpha) { # test passed
            message("Chi2 test PASSed with p-value=",signif(pval))
            break
          }
        
          else {
            if (i<1.5) 
              warning("Chi2 test failed once with p-value=",signif(pval)," --> Try again")
            else
              stop("chi2 test FAILED!  p-value=",signif(pval), call.=FALSE)
          }
        }

} ## --- end of rstream.check.chi2() ---

## ..........................................................................

rstream.check.chi2.run <- function (type, ...) {
        ##  Run a chi^2 test and evaluate p-value.
        ##  
        ##  type   ...  type (class) of stream 
        ##  ...    ...  optional args for new(...)
        ##
        
	## -- Create a stream 
	s <- new (type, ...)
        print(s)

        ## -- Draw sample 
	u <- rstream.sample(s,chi2.samplesize)
	if (length(u) != chi2.samplesize)
		stop("chi2: sample.rstream failed")

        ## -- make histogram of with equalsized bins (classified data)
        nbins <- as.integer(sqrt(chi2.samplesize))
        breaks <- (0:nbins)/nbins
        h <- hist(u,plot=F,breaks=breaks)$count
        ## -- run unur.chiq.test
        pval <- chisq.test(h)$p.value

        ## -- return p-value
        pval

} ## --- end of rstream.check.chi2.run() ---


#############################################################################
##                                                                          #
##  Auxiliary routines                                                      #
##                                                                          #
#############################################################################

## Test whether there is an error -------------------------------------------

iserror <- function (expr) { is(try(expr), "try-error") }

## -- End -------------------------------------------------------------------
