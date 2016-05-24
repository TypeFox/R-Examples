
stopifnot(require(RDieHarder, quiet=TRUE))

dh <- dieharder("rand", "diehard_runs", seed=12345)

print(summary(dh))
