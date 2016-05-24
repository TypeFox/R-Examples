#!/usr/bin/r

library(RVowpalWabbit)

## change to 'test' directory of package
setwd( system.file("test", package="RVowpalWabbit") )


# Test 6: run predictions on Test 4 model
# Pretending the labels aren't there
# {VW} -t -i models/0002.model -d train-sets/0002.dat -p 0002b.predict
test6 <- c("-t", "-i", "models/0002.model",
           "-d", "train-sets/0002.dat",
           "-p", "/tmp/0002b.predict")

res <- vw(test6, quiet=FALSE)
print(res)
