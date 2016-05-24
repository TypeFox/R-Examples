source("local.R")

if (test) { 

library("BRugs")

## Prepare the example files in a temporary directory
exfiles <- dir(options()$OpenBUGSExamples, pattern="^Rats.*txt$", full.names=TRUE)
ok <- file.copy(exfiles, tempdir())
if(!all(ok)) 
    stop("Some files could not be copied from OpenBUGS examples to the temporary directory")

exfiles <- dir(options()$OpenBUGSExamples, pattern="^Beetles.*txt$", full.names=TRUE)
ok <- file.copy(exfiles, tempdir())
if(!all(ok)) 
    stop("Some files could not be copied from OpenBUGS examples to the temporary directory")


BRugsFit(data = "Ratsdata.txt", inits = "Ratsinits.txt",
    para = c("alpha", "beta"), modelFile = "Ratsmodel.txt",
    numChains = 1,
    working.directory = tempdir())

setwd(tempdir())
modelCheck("Ratsmodel.txt")
modelData("Ratsdata.txt")
modelCompile(numChains=2)
modelInits(rep("Ratsinits.txt", 2))
modelUpdate(1000)
samplesSet(c("alpha0", "alpha"))
modelUpdate(1000)
samplesStats("*")

### Four different ways of supplying the data 

beetles <- list(x = c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839),
                n = c(59, 60, 62, 56, 63, 59, 62, 60),
                r = c(6, 13, 18, 28, 52, 53, 61, 60), N = 8)

BRugsFit(data = "Beetlesdata.txt", inits = "Beetlesinits.txt",
         para = c("alpha", "beta", "rhat"), modelFile = "Beetlesmodel.txt",
         numChains = 1,
         working.directory = tempdir())

BRugsFit(data = beetles, inits = "Beetlesinits.txt",
         para = c("alpha", "beta", "rhat"), modelFile = "Beetlesmodel.txt",
         numChains = 1,
         working.directory = tempdir())

with(beetles, 
     BRugsFit(data = list("x", "n", "r", "N"), inits = "Beetlesinits.txt",
              para=c("alpha", "beta", "rhat"), modelFile = "Beetlesmodel.txt",
              numChains = 1, working.directory = tempdir())
     )

with(beetles, 
     BRugsFit(data = c("x", "n", "r", "N"), inits = "Beetlesinits.txt",
              para=c("alpha", "beta", "rhat"), modelFile = "Beetlesmodel.txt",
              numChains = 1, working.directory = tempdir())
     )

}
