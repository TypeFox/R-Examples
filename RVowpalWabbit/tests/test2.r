#!/usr/bin/r

library(RVowpalWabbit)

## change to 'test' directory of package
setwd( system.file("test", package="RVowpalWabbit") )


# t 2: checking predictions as well
# {VW} -t train-sets/0001.dat -i models/0001.model -p 001.predict.tmp
test2 <- c("-t", "train-sets/0001.dat",
           "-i", "models/0001.model",
           "-p", "/tmp/0001.predict.tmp")

res <- vw(test2, quiet=FALSE)
print(res)
