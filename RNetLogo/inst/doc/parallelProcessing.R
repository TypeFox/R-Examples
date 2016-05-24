### R code from vignette source 'parallelProcessing.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: simpleFunction
###################################################
testfun1 <- function(x) {
  return(x*x)
}


###################################################
### code chunk number 3: applySimpleFunction
###################################################
my.v1 <- 1:10
print(my.v1)

my.v1.quad <- sapply(my.v1, testfun1)
print(my.v1.quad)


###################################################
### code chunk number 4: makeCluster1
###################################################
# load the parallel package
library(parallel)

# detect the number of cores available
processors <- detectCores()


###################################################
### code chunk number 5: parallelProcessing.Rnw:133-136
###################################################
if (processors > 2) {
  processors <- 2
}


###################################################
### code chunk number 6: parallelProcessing.Rnw:138-140
###################################################
# create a cluster
cl <- makeCluster(processors)


###################################################
### code chunk number 7: runCluster1
###################################################
# call parallel sapply
my.v1.quad.par <- parSapply(cl, my.v1, testfun1)
print(my.v1.quad.par)

# stop cluster
stopCluster(cl)


###################################################
### code chunk number 8: defineFunctions1
###################################################
# the initialization function
prepro <- function(dummy, gui, nl.path, model.path) {
  library(RNetLogo)
  NLStart(nl.path, gui=gui)
  NLLoadModel(model.path)
}


# the simulation function
simfun <- function(x) {
  NLCommand("print ",x)
  NLCommand("set density", x)
  NLCommand("setup")
  NLCommand("go")
  NLCommand("print count turtles")
  ret <- data.frame(x, NLReport("count turtles"))
  names(ret) <- c("x","turtles")
  return(ret)
}


# the quit function
postpro <- function(x) {
  NLQuit()
}


###################################################
### code chunk number 9: defineFunctions1XX
###################################################
# the initialization function
prepro <- function(dummy, gui, nl.path, model.path) {
  library(RNetLogo)
  NLStart(nl.path, gui=FALSE)
  NLLoadModel(model.path)
}


###################################################
### code chunk number 10: initNetLogo1 (eval = FALSE)
###################################################
## # load the parallel package, if not already done
## require(parallel)
## 
## # detect the number of cores available
## processors <- detectCores()
## 
## # create cluster
## cl <- makeCluster(processors)
## 
## # set variables for the start up process
## # adapt path appropriate (or set an environment variable NETLOGO_PATH)
## gui <- TRUE
## nl.path <- Sys.getenv("NETLOGO_PATH", "C:/Program Files/NetLogo 5.3/app")
## model.path <- "models/Sample Models/Earth Science/Fire.nlogo"
## 
## # load NetLogo in each processor/core
## invisible(parLapply(cl, 1:processors, prepro, gui=gui, 
##                     nl.path=nl.path, model.path=model.path))


###################################################
### code chunk number 11: simNetLogo1
###################################################
# create a vector with 20 density values
density <- 1:20
print(density)


###################################################
### code chunk number 12: simNetLogo1b (eval = FALSE)
###################################################
## # run a simulation for each density value
## # by calling parallel sapply
## result.par <- parSapply(cl, density, simfun)


###################################################
### code chunk number 13: parallelProcessing.Rnw:241-242 (eval = FALSE)
###################################################
## save.image(file="parallelprocessData1.RData")


###################################################
### code chunk number 14: parallelProcessing.Rnw:244-245
###################################################
load("parallelprocessData1.RData")


###################################################
### code chunk number 15: parallelProcessing.Rnw:247-248
###################################################
print(data.frame(t(result.par)))


###################################################
### code chunk number 16: quitNetLogo1 (eval = FALSE)
###################################################
## # Quit NetLogo in each processor/core
## invisible(parLapply(cl, 1:processors, postpro))
## 
## # stop cluster
## stopCluster(cl)


###################################################
### code chunk number 17: simNetLogo1 (eval = FALSE)
###################################################
## # run in headless mode
## gui <- FALSE
## 
## # create cluster
## cl <- makeCluster(processors)
## 
## # load NetLogo in each processor/core
## invisible(parLapply(cl, 1:processors, prepro, gui=gui, 
##                     nl.path=nl.path, model.path=model.path))


###################################################
### code chunk number 18: parallelProcessing.Rnw:280-283
###################################################
# create a vector with 20 density values
density <- 1:20
print(density)


###################################################
### code chunk number 19: parallelProcessing.Rnw:285-288 (eval = FALSE)
###################################################
## # run a simulation for each density value
## # by calling parallel sapply
## result.par <- parSapply(cl, density, simfun)


###################################################
### code chunk number 20: parallelProcessing.Rnw:290-291 (eval = FALSE)
###################################################
## save.image(file="parallelprocessData2.RData")


###################################################
### code chunk number 21: parallelProcessing.Rnw:293-294
###################################################
load("parallelprocessData2.RData")


###################################################
### code chunk number 22: parallelProcessing.Rnw:296-297
###################################################
print(data.frame(t(result.par)))


###################################################
### code chunk number 23: parallelProcessing.Rnw:299-304 (eval = FALSE)
###################################################
## # Quit NetLogo in each processor/core
## invisible(parLapply(cl, 1:processors, postpro))
## 
## # stop cluster
## stopCluster(cl)


