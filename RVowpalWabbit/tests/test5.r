#!/usr/bin/r

library(RVowpalWabbit)

## change to 'test' directory of package
setwd( system.file("test", package="RVowpalWabbit") )

# Test 5: add -q .., adaptive, and more (same input, different outputs)
# {VW} --initial_t 1 --power_t 0.5 --adaptive -q Tf -q ff -f models/0002a.model train-sets/0002.dat
test5 <- c("--initial_t", "1",
           "--power_t", "0.5",
           "--adaptive",
           "-q", "Tf",
           "-q", "ff",
           "-f", "models/0002a.model",
           "train-sets/0002.dat")

res <- vw(test5, quiet=FALSE)
print(res)
