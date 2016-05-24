#############################################################################
##                                                                         ##
##   Tests for class: rstream.runif                                        ##
##                                                                         ##
#############################################################################

## Test Parameters ----------------------------------------------------------

## samplesize <- 1e4
## streamname <- "testname"


## Load test routines -------------------------------------------------------

source ("test_routines.R")


## Check: API ---------------------------------------------------------------

## Run general tests
rstream.check.API("rstream.runif", kind="current", t_incp=FALSE)
rstream.check.API("rstream.runif", kind="default", t_incp=FALSE)
rstream.check.API("rstream.runif", kind="Wichmann-Hill", t_incp=FALSE)
rstream.check.API("rstream.runif", kind="Marsaglia-Multicarry", t_incp=FALSE)
rstream.check.API("rstream.runif", kind="Super-Duper", t_incp=FALSE)
rstream.check.API("rstream.runif", kind="Mersenne-Twister", t_incp=FALSE)
rstream.check.API("rstream.runif", kind="Knuth-TAOCP", t_incp=FALSE)
rstream.check.API("rstream.runif", kind="Knuth-TAOCP-2002", t_incp=FALSE)

## user-defined RNG must not work
if( ! iserror( rstream.check.API("rstream.runif", kind="user", t_incp=FALSE) ) )
	stop("\"user-defined\" must not work")

## Test interface rstream <-> R RNG
rstream.check.setRNG("rstream.runif", kind="current")
rstream.check.setRNG("rstream.runif", kind="default")
rstream.check.setRNG("rstream.runif", kind="Wichmann-Hill")
rstream.check.setRNG("rstream.runif", kind="Marsaglia-Multicarry")


## Special tests -----------------------------------------------------------

## run with all optional arguments
new("rstream.runif",seed=12345, antithetic=TRUE)

## when using kind="current" (which is the default)
## runif() and rstream.sample should return the same values (until rstream.reset)
s <- new("rstream.runif")
x <- rstream.sample(s,samplesize)
y <- runif(samplesize)
if (!identical(all.equal(x, y), TRUE))
	stop("\"current\" failed:1")

## rstream.sample() must not change .Random.seed
RNGkind("Knuth-TAOCP")
save.seed <- .Random.seed
x <- rstream.sample(s,samplesize)
if (!identical(all.equal(save.seed, .Random.seed), TRUE))
	stop("seed not restored after calling rstream.sample()")


## End ----------------------------------------------------------------------
