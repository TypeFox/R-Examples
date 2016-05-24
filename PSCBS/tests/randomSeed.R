library("PSCBS")

message("*** randomSeed() - setup ...")
ovars <- ls(envir=globalenv())
genv <- globalenv()
RNGkind("Mersenne-Twister")
if (exists(".Random.seed", envir=genv, inherits=FALSE))
  rm(list=".Random.seed", envir=genv, inherits=FALSE)
seed0 <- genv$.Random.seed
stopifnot(is.null(seed0))
okind0 <- RNGkind()[1L]

sample1 <- function() { sample(0:9, size=1L) }
message("*** randomSeed() - setup ... done")


message("*** randomSeed('get') ...")
## Get random seed
seed <- randomSeed("get")
stopifnot(identical(seed, seed0))

## Repeat after new sample
y1 <- sample1()
message(sprintf("Random number: %d", y1))
seed1 <- randomSeed("get")
stopifnot(!identical(seed1, seed0))
message("*** randomSeed('get') ... done")


message("*** randomSeed('set', 42L) ...")
randomSeed("set", seed=42L)
seed2 <- randomSeed("get")
stopifnot(!identical(seed2, seed1))

y2 <- sample1()
message(sprintf("Random number: %d (with random seed = 42L)", y2))

## Reset to previous state
randomSeed("reset")
seed3 <- randomSeed("get")
stopifnot(identical(seed3, seed1))
stopifnot(identical(RNGkind()[1L], okind0),
          identical(randomSeed("get"), seed1))
message("*** randomSeed('set', 42L) ... done")


message("*** randomSeed('set', NULL) ...")
randomSeed("set", seed=NULL)
seed4 <- randomSeed("get")
stopifnot(is.null(seed4))

y3 <- sample1()
message(sprintf("Random number: %d", y3))

message("*** randomSeed('set', NULL) ... done")


message("*** randomSeed('set', 42L) again ...")
seed5 <- randomSeed("get")
randomSeed("set", seed=42L)
y4 <- sample1()
message(sprintf("Random number: %d (with random seed = 42L)", y4))
stopifnot(identical(y4, y2))

randomSeed("reset")
stopifnot(identical(RNGkind()[1L], okind0),
          identical(randomSeed("get"), seed5))
message("*** randomSeed('set', 42L) again ... done")



## L'Ecuyer-CMRG: Random number generation for parallel processing
message("*** randomSeed(): L'Ecuyer-CMRG stream ...")

okind <- RNGkind()[1L]
stopifnot(identical(okind, okind0))

randomSeed("set", seed=NULL)
oseed <- randomSeed("get")
stopifnot(is.null(oseed))

randomSeed("set", seed=42L, kind="L'Ecuyer-CMRG")
oseed2 <- randomSeed("reset")
str(oseed2)
stopifnot(identical(oseed2, oseed))
stopifnot(identical(RNGkind()[1L], okind),
          identical(randomSeed("get"), oseed))

randomSeed("set", seed=42L, kind="L'Ecuyer-CMRG")
seed0 <- randomSeed("get")
seeds0 <- lapply(1:10, FUN=function(i) randomSeed("advance"))
oseed2 <- randomSeed("reset")
stopifnot(identical(oseed2, oseed))
stopifnot(identical(RNGkind()[1L], okind),
          identical(randomSeed("get"), oseed))


## Assert reproducible .Random.seed stream
randomSeed("set", seed=42L, kind="L'Ecuyer-CMRG")
seed1 <- randomSeed("get")
seeds1 <- lapply(1:10, FUN=function(i) randomSeed("advance"))
stopifnot(identical(seed1, seed0))
stopifnot(identical(seeds1, seeds0))

randomSeed("reset")
stopifnot(identical(RNGkind()[1L], okind),
          identical(randomSeed("get"), oseed))

randomSeed("set", seed=42L, kind="L'Ecuyer-CMRG")
seeds2 <- randomSeed("advance", n=10L)
stopifnot(identical(seeds2, seeds0))

randomSeed("reset")
stopifnot(identical(RNGkind()[1L], okind),
          identical(randomSeed("get"), oseed))

randomSeed("set", seed=seeds2[[1]], kind="L'Ecuyer-CMRG")
randomSeed("reset")
stopifnot(identical(RNGkind()[1L], okind),
          identical(randomSeed("get"), oseed))

randomSeed("set", seed=42L, kind="L'Ecuyer-CMRG")
y0 <- sapply(1:10, FUN=function(ii) {
  randomSeed("advance")
  sample1()
})
print(y0)
randomSeed("reset")

randomSeed("set", seed=42L, kind="L'Ecuyer-CMRG")
y1 <- sapply(1:10, FUN=function(ii) {
  randomSeed("advance")
  sample1()
})
print(y1)
stopifnot(identical(y1, y0))
randomSeed("reset")

stopifnot(identical(RNGkind()[1L], okind))

message("*** randomSeed(): L'Ecuyer-CMRG stream ... done")


## Cleanup
message("*** randomSeed() - cleanup ...")
genv <- globalenv()
RNGkind("Mersenne-Twister")
if (exists(".Random.seed", envir=genv, inherits=FALSE))
  rm(list=".Random.seed", envir=genv, inherits=FALSE)
rm(list=ovars, envir=globalenv())
message("*** randomSeed() - cleanup ... done")
