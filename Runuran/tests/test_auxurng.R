
#############################################################################
##                                                                         ##
##   Tests for class 'unuran'                                              ##
##                                                                         ##
#############################################################################

## --- Test Parameters ------------------------------------------------------

## size of sample for test
samplesize <- 1.e5

## --- Load library ---------------------------------------------------------

library(Runuran)

## --------------------------------------------------------------------------

## Test whether there is an error
iserror <- function (expr) { is(try(expr), "try-error") }

## --------------------------------------------------------------------------

if (! iserror( set.aux.seed()))   stop("missing seed not detected")
if (! iserror( set.aux.seed(0)))  stop("seed=0 not detected")
if (! iserror( set.aux.seed(-1))) stop("seed < 0 not detected")

## --------------------------------------------------------------------------

gen <- unuran.new(udnorm(), "srou");
if (! is.na(use.aux.urng(gen))) stop("use.aux.urng() must return NA")
rm(gen)

## --------------------------------------------------------------------------

gen1 <- unuran.new(udnorm(), "tdr; cpoints=2; max_sqhratio=0.5; usedars=on")
gen2 <- unuran.new(udexp(), "tdr; cpoints=2; max_sqhratio=0.5; usedars=on")

set.seed(12345)
x1 <- ur(gen1,samplesize)
set.seed(12345)
x2 <- ur(gen2,samplesize)
if (abs(cor(x1,x2)) > 0.5) stop("correlation too high: aux urng must not be used")

if (use.aux.urng(gen1)) stop("use.aux.urng() must return FALSE")

use.aux.urng(gen1) <- TRUE
use.aux.urng(gen2) <- TRUE

set.seed(12345)
x1 <- ur(gen1,samplesize)
set.seed(12345)
x2 <- ur(gen2,samplesize)
if (abs(cor(x1,x2)) < 0.5) stop("correlation too small: aux urng should be used")

if (! use.aux.urng(gen1)) stop("use.aux.urng() must return TRUE")

use.aux.urng(gen1) <- FALSE
use.aux.urng(gen2) <- FALSE

set.seed(12345)
x1 <- ur(gen1,samplesize)
set.seed(12345)
x2 <- ur(gen2,samplesize)
if (abs(cor(x1,x2)) > 0.5) stop("correlation too high: aux urng must not be used")

if (use.aux.urng(gen1)) stop("use.aux.urng() must return FALSE")

rm(gen1,gen2)

## --------------------------------------------------------------------------

gen <- unuran.new(udnorm(), "tdr; cpoints=2; max_sqhratio=0.5; usedars=on")
use.aux.urng(gen) <- TRUE

set.seed(12345)
x1 <- ur(gen,samplesize)
set.seed(12345)
x2 <- ur(gen,samplesize)

if (isTRUE(all.equal(x1,x2))) stop("strings should differ: aux urng not working properly")

set.seed(12345)
set.aux.seed(999)
x1 <- ur(gen,samplesize)
set.seed(12345)
set.aux.seed(999)
x2 <- ur(gen,samplesize)

if (!isTRUE(all.equal(x1,x2))) stop("strings must not differ: aux urng not working properly")

rm(gen)

## --- End ------------------------------------------------------------------

detach("package:Runuran",unload = TRUE)

## --------------------------------------------------------------------------
