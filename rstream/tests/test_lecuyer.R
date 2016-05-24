#############################################################################
##                                                                         ##
##   Tests for class: rstream.lecuyer                                      ##
##                                                                         ##
#############################################################################

## Test Parameters ----------------------------------------------------------

## samplesize <- 1e4
## streamname <- "testname"


## Load test routines -------------------------------------------------------

source ("test_routines.R")


## Check: API ---------------------------------------------------------------

## Run general tests
rstream.check.API("rstream.lecuyer")

## Test interface rstream <-> R RNG
rstream.check.setRNG("rstream.lecuyer")


## Check: Goodness-of-fit ---------------------------------------------------

## Chi^2 goodness-of-fit test
rstream.check.chi2("rstream.lecuyer")
rstream.check.chi2("rstream.lecuyer", antithetic=TRUE)
rstream.check.chi2("rstream.lecuyer", incprecision=TRUE)


## Special tests ------------------------------------------------------------

## it must not be possible to reset the seed without flag 'force.seed'
if( ! iserror( new("rstream.lecuyer", seed=rep(12345,6)) ) )
	stop("force.seed = TRUE required")

## one must give at least 6 numbers for the seed
if( ! iserror( new("rstream.lecuyer", seed=rep(12345,5), force.seed=TRUE) ) )
	stop("seed must have length 6")

## run with all optional arguments
s <- new("rstream.lecuyer", seed=rep(12345,6), force.seed=TRUE, antithetic=TRUE, incprecision=TRUE)

## check generator
check_gen <-
        function(s, mean, var) {
                x <- rstream.sample(s,100)
                m <- mean(x)
                v <- 99/100 * var(x)
                if( abs(mean-m) > 1.e-9) stop ("wrong mean")
                if( abs(var-v) > 1.e-9) stop ("wrong variance")
}
	
s1 <- new("rstream.lecuyer", seed=rep(12345,6), force.seed=TRUE)
s2 <- new("rstream.lecuyer")
s3 <- new("rstream.lecuyer")

check_gen(s=s1, mean=0.51783564702603, var=0.08251698866776)
check_gen(s=s2, mean=0.52394642023622, var=0.09989286188901)
check_gen(s=s3, mean=0.54064030708121, var=0.09043368375419)

check_gen(s=s2, mean=0.45591975817022, var=0.07421246146254)
rstream.reset(s2)
check_gen(s=s2, mean=0.52394642023622, var=0.09989286188901)

rstream.antithetic(s2) <- TRUE
check_gen(s=s2, mean=0.54408024182978, var=0.07421246146254)

rstream.antithetic(s2) <- FALSE
rstream.incprecision(s2) <- TRUE
check_gen(s=s2, mean=0.57325868991301, var=0.07898987147185)

rstream.antithetic(s2) <- TRUE
check_gen(s=s2, mean=0.56487612596530, var=0.08846386734202)


## End ----------------------------------------------------------------------
