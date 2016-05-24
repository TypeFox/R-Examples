#!/usr/bin/r

library(RVowpalWabbit)

## change to 'test' directory of package
setwd( system.file("test", package="RVowpalWabbit") )


# Test 7: using -q and multiple threads
# {VW} --adaptive -q ff -f models/0002c.model train-sets/0002.dat
test7 <- c("--adaptive", "-q", "ff",
           "-f", "models/0002c.model",
           "train-sets/0002.dat")

res <- vw(test7, quiet=FALSE)
print(res)
