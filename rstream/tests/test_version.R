## --------------------------------------------------------------------------
## Test loading of saved rstream object
##

library(rstream)

## -- Bug fixes in 1.3 ------------------------------------------------------

## sample size for test
samplesize <- 10

## -- rstream.mrg32k3a ......................................................

## create stream
s <- new("rstream.mrg32k3a")

## make a working copy
sc <- rstream.clone(s)

## get sequence
x <- rstream.sample(sc, samplesize)

## set RNG
rstream.RNG(s)

## compare sequence
y <- runif(samplesize)
if (!identical(all.equal(x, y), TRUE))
  stop("set version failed for rstream.mrg32k3a: 1")

## try again
rstream.RNG(s)
y <- runif(samplesize)
if (!identical(all.equal(x, y), TRUE))
  stop("set version failed for rstream.mrg32k3a: 2")

## now switch version
rstream.version("1.2")
rstream.version()

rstream.RNG(s)
y <- runif(samplesize)
if (identical(all.equal(x, y), TRUE))
  stop("set version failed for rstream.mrg32k3a: 3")

if (! identical(all.equal(x[2:samplesize], y[1:(samplesize-1)]), TRUE))
  stop("set version failed for rstream.mrg32k3a: 4")

## lastest version
rstream.version("default")
rstream.version()

rstream.RNG(s)
y <- runif(samplesize)
if (!identical(all.equal(x, y), TRUE))
  stop("set version failed for rstream.mrg32k3a: 5")


## -- rstream.lecuyer ......................................................

## Remark: this class is deprecated

## create stream
s <- new("rstream.lecuyer")

## make a working copy
sc <- rstream.clone(s)

## get sequence
x <- rstream.sample(sc, samplesize)

## set RNG
rstream.RNG(s)

## compare sequence
y <- runif(samplesize)
if (!identical(all.equal(x, y), TRUE))
  stop("set version failed for rstream.lecuyer: 1")

## try again
rstream.RNG(s)
y <- runif(samplesize)
if (!identical(all.equal(x, y), TRUE))
  stop("set version failed for rstream.lecuyer: 2")

## now switch version
rstream.version("1.2")
rstream.version()

rstream.RNG(s)
y <- runif(samplesize)
if (identical(all.equal(x, y), TRUE))
  stop("set version failed for rstream.lecuyer: 3")

if (! identical(all.equal(x[2:samplesize], y[1:(samplesize-1)]), TRUE))
  stop("set version failed for rstream.lecuyer: 4")

## lastest version
rstream.version("default")
rstream.version()

rstream.RNG(s)
y <- runif(samplesize)
if (!identical(all.equal(x, y), TRUE))
  stop("set version failed for rstream.lecuyer: 5")


## -- End -------------------------------------------------------------------
