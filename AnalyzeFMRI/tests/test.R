###################################
## Tests for AnalyzeFMRI package ##
###################################

library(AnalyzeFMRI)

## Test that reading of examples works

a1 <- f.read.analyze.volume(system.file("example.img", package="AnalyzeFMRI"))

a1[30:40, 30:40, 10, 1]


## Test that writing of examples works

a2 <- array(1:1000, dim = c(10, 10, 10, 1))
f.write.analyze(a2, file = "test.array", size = "float")

a3 <- f.read.analyze.volume(file = "test.array.img")

sum(a2 != a3) ## should return 0

