source("local.R")

if (test) {

library("BRugs")
exfiles <- dir(options()$OpenBUGSExamples, pattern="^Rats.*txt$", full.names=TRUE)
ok <- file.copy(exfiles, tempdir())
if(!all(ok)) 
    stop("Some files could not be copied from OpenBUGS examples to the temporary directory")
setwd(tempdir())
## .onLoad(lib=.libPaths()[1], pkg="BRugs") # if developing without using namespace

### TEST ALL USER-LEVEL FUNCTIONS USING RATS EXAMPLE

### QUESTIONS
### any need to worry about executable files: .dlls in source package ?
### why not download all the openbugs stuff?  would need a static webpage
### Do we still support S-Plus?
### How do we restrict it to Linux and Windows, and stop it from trying to compile on e.g. Mac/Solaris i386?
### FIXME
# setwd("/home/chris/openbugs/gappy dir") # works
# setwd("/home/chris/openbugs/\'quotey dir\'") # doesn't work
# setwd("/home/chris/openbugs/\"quotey dir\"") # doesn't work

### Basic model setup and stats
modelCheck("Ratsmodel.txt")
modelData("Ratsdata.txt")
modelCompile(numChains=2)
modelInits("Ratsinits.txt", 1)
modelInits("Ratsinits1.txt", 2)
#modelInits("Ratsinits_noran.txt", 1)
#modelInits("Ratsinits1_noran.txt", 2)
#modelGenInits()
modelSetRN(4)
modelUpdate(1000)
samplesSet(c("alpha0", "alpha"))
summarySet(c("alpha0","alpha"))
ranksSet("alpha")
dicSet()
modelUpdate(1000)
stopifnot(samplesMonitors("*")[31]=="alpha0")
stopifnot(isTRUE(all.equal(samplesStats("*")["alpha0","mean"], 106.6)))
stopifnot(isTRUE(all.equal(dicStats()["total","DIC"],1021)))
stopifnot(isTRUE(all.equal(summaryStats("alpha0")["alpha0","mean"], 106.6)))
stopifnot(isTRUE(all.equal(as.numeric(ranksStats("alpha")[1,]), c(10,12,18))))

### Utilities
samplesSetBeg(1200); stopifnot(samplesGetBeg()==1200)
samplesSetEnd(1500); stopifnot(samplesGetEnd()==1500)
samplesSetFirstChain(2); stopifnot(samplesGetFirstChain()==2)
samplesSetLastChain(2); stopifnot(samplesGetLastChain()==2)
samplesSetThin(2); stopifnot(samplesGetThin()==2)
samplesSetBeg(1000)
samplesSetEnd(1000000)
samplesSetFirstChain(1)
samplesSetLastChain(2)
samplesSetThin(1)
stopifnot(samplesSize("alpha0")==2*1000)
stopifnot(length(samplesSample("alpha0"))==2*1000)
stopifnot(all(dim(infoNodeValues("alpha"))==c(30,2)))
stopifnot(BRugs:::dimensions("alpha")==1)
stopifnot(BRugs:::dimensions("alpha0")==0)
stopifnot(modelIteration()==2000)
stopifnot(all.equal(sort(modelNames()),sort(c("N", "T", "Y", "alpha", "alpha.c", "alpha.tau", "alpha0", "beta","beta.c", "beta.tau", "deviance", "mu", "sigma", "tau.c", "x","xbar"))))
modelPrecision(8)
stopifnot(nchar(as.character(samplesStats("alpha0")$mean))==9)
modelPrecision(4)
x <- infoNodeValues("alpha")
setValues("alpha", x-1)
stopifnot(isTRUE(all.equal(infoNodeValues("alpha"), x-1, tol=1e-06)))
stopifnot(modelAdaptivePhase()==-1)

### Plots
stopifnot(all(dim(samplesHistory("alpha0", plot=interactive())[[1]])==c(1000,2)))
stopifnot(isTRUE(all.equal(samplesAutoC("alpha0", 1, plot=interactive())$alpha0$acf[1], 1)))
stopifnot(all(dim(samplesBgr("alpha0", plot=interactive())$alpha0)==c(50,4)))
stopifnot(all.equal(samplesBgr("alpha0", plot=FALSE)$alpha0$pooled[1], 0.411, tol=0.1))
# stopifnot(isTRUE(all.equal(samplesCorrel("alpha[1]", "alpha[2]")[1,1], 0.8924, tol=0.1)))
stopifnot(all(dim(samplesDensity("alpha", plot=interactive(), ask=FALSE))==c(7,30)))
stopifnot(all(dim(plotHistory("alpha0", plot=interactive())[[1]])==c(1000,2)))
stopifnot(isTRUE(all.equal(plotAutoC("alpha0", 1, plot=interactive())$acf[1], 1)))
stopifnot(all(dim(plotBgr("alpha0", plot=interactive()))==c(50,4)))
stopifnot(length(plotDensity("alpha0",plot=FALSE))==7)

### Clearing
samplesClear("alpha")
stopifnot(samplesMonitors("*")=="alpha0")
stopifnot(samplesSize("alpha")==0)
dicClear()
summaryClear("*")

## External and file access
samplesCoda("*" ,stem="Rats")
modelSaveState(stem="Rats")
require("coda")
rats.coda <- buildMCMC("alpha0")
stopifnot(length(rats.coda[[1]])==1000)
modelCheck("Ratsmodel.txt")
dat <- dget("Ratsdata.txt")
modelData(bugsData(dat))
modelCompile(numChains=2)
inits <- dget("Ratsinits.txt")
modelInits(bugsInits(list(inits,inits), numChains=2))
# writeModel() # tested in help page example

## Manuals
if (interactive()) help.BRugs() # todo remove link from doc
if (interactive()) help.WinBUGS()

## Gen inits
modelCheck("Ratsmodel.txt")
modelData("Ratsdata.txt")
modelCompile(numChains=2)
inits <- dget("Ratsinits.txt"); inits[c("alpha","beta")] <- NULL
modelInits(bugsInits(list(inits,inits), numChains=2))
modelGenInits()

## Functions to tune updaters.  Just make sure their interfaces are correct
#modelSetAP()
#modelSetIts()
#modelSetOR()
#modelEnable()
#modelDisable()

## Info functions
stopifnot(infoNodeMethods("alpha")[,"Type"] == "UpdaterNormal.StdUpdater")
stopifnot(infoNodeTypes("alpha")[1,"Type"]=="GraphNormal.Node")
stopifnot(infoUpdatersbyName()["alpha.c","Type"]=="conjugate normal updater")
stopifnot(infoUpdatersbyDepth()["alpha.c","Type"]=="conjugate normal updater")
mem <- infoMemory()
stopifnot(is.numeric(mem) && mem > 0)
stopifnot(any(infoModules()[,"Module"] == "Kernel"))
}
