#!/usr/bin/r

library(RVowpalWabbit)

## change to 'test' directory of package
setwd( system.file("test", package="RVowpalWabbit") )


# Test 1:
# {VW} -b 17 -l 20 --initial_t 128000 --power_t 1 -d train-sets/0001.dat -f models/0001.model -c --passes 2 --compressed --ngram 3 --skips 1
test1 <- c("-b", "17",
           "-l", "20",
           "--initial_t", "128000",
           "--power_t", "1",
           "-d", "train-sets/0001.dat",
           "-f", "models/0001.model",
           "-c",
           "--passes", "2",
           "--compressed",
           "--ngram", "3",
           "--skips", "1")

res <- vw(test1, quiet=FALSE)
print(res)
