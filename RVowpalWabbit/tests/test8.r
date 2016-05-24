#!/usr/bin/r

library(RVowpalWabbit)

## change to 'test' directory of package
setwd( system.file("test", package="RVowpalWabbit") )


## Test 8: predicts on test 7 model
## {VW} -t -i models/0002c.model -d train-sets/0002.dat -p 0002c.predict
test8 <- c("-t",
           "-i", "models/0002c.model",
           "-d", "train-sets/0002.dat",
           "-p", "/tmp/0002c.predict")

res <- vw(test8, quiet=FALSE)
print(res)
