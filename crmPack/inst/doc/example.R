### R code from vignette source 'example.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: load
###################################################
## load package already here to be able to reference the package version!!
library(crmPack)
crmPackVersion <- as.character(packageVersion("crmPack"))


###################################################
### code chunk number 2: setup
###################################################
options(continue="  ")                  # use two blanks instead of "+" for
                                        # continued lines (for easy copying)


###################################################
### code chunk number 3: load
###################################################
library(crmPack)


###################################################
### code chunk number 4: webinterface (eval = FALSE)
###################################################
## crmPackHelp()


###################################################
### code chunk number 5: emptydata
###################################################
emptydata <- Data(doseGrid=
                 c(0.1, 0.5, 1.5, 3, 6,
                   seq(from=10, to=80, by=2)))


###################################################
### code chunk number 6: data
###################################################
data <- Data(x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
             y=c(0, 0, 0, 0, 0, 0, 1, 0),
             cohort=c(0, 1, 2, 3, 4, 5, 5, 5),
             doseGrid=
                 c(0.1, 0.5, 1.5, 3, 6,
                   seq(from=10, to=80, by=2)))


###################################################
### code chunk number 7: ids
###################################################
data@ID


###################################################
### code chunk number 8: plotdata
###################################################
print(plot(data))


###################################################
### code chunk number 9: model-setup
###################################################
model <- LogisticLogNormal(mean=c(-0.85, 1),
                           cov=
                               matrix(c(1, -0.5, -0.5, 1),
                                      nrow=2),
                           refDose=56)


###################################################
### code chunk number 10: class
###################################################
class(model)


###################################################
### code chunk number 11: str
###################################################
str(model)


###################################################
### code chunk number 12: dose
###################################################
model@dose


###################################################
### code chunk number 13: emptydata
###################################################
emptydata <- Data(doseGrid=
                   seq(from=25, to=300, by=25))
data1 <- emptydata                 


###################################################
### code chunk number 14: LogisticIndepBeta
###################################################
DLTmodel<-LogisticIndepBeta(binDLE=c(1.05,1.8),DLEweights=c(3,3),
                            DLEdose=c(25,300),data=data1)


###################################################
### code chunk number 15: strLogisticIndepBeta
###################################################
str(DLTmodel)


###################################################
### code chunk number 16: phi1
###################################################
DLTmodel@phi1


###################################################
### code chunk number 17: min-inf
###################################################
coarseGrid <- c(0.1, 10, 30, 60, 100)
minInfModel <- MinimalInformative(dosegrid = coarseGrid,
                                  refDose=50,
                                  threshmin=0.2,
                                  threshmax=0.3,
                                  control=
                                 list(threshold.stop=0.03,
                                       maxit=200),
                                  seed=432)


###################################################
### code chunk number 18: min-inf-res
###################################################
matplot(x=coarseGrid,
        y=minInfModel$required,
        type="b", pch=19, col="blue", lty=1,
        xlab="dose",
        ylab="prior probability of DLT")
matlines(x=coarseGrid,
         y=minInfModel$quantiles,
         type="b", pch=19, col="red", lty=1)
legend("right",
       legend=c("quantiles", "approximation"),
       col=c("blue", "red"),
       lty=1,
       bty="n")


###################################################
### code chunk number 19: min-inf-dist
###################################################
minInfModel$distance


###################################################
### code chunk number 20: min-inf-model
###################################################
str(minInfModel$model)


###################################################
### code chunk number 21: min-inf-model-extract
###################################################
myModel <- minInfModel$model


###################################################
### code chunk number 22: mcmc-opts
###################################################
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=2000)


###################################################
### code chunk number 23: mcmc-sampling
###################################################
set.seed(94)
samples <- mcmc(data, model, options)


###################################################
### code chunk number 24: mcmc-extract
###################################################
## look at the structure of the samples object:
str(samples)
## now extract the alpha0 samples (intercept of the regression model)
alpha0samples <- get(samples, "alpha0")


###################################################
### code chunk number 25: ggmcmc
###################################################
library(ggmcmc)
print(ggs_traceplot(alpha0samples))


###################################################
### code chunk number 26: ggmcmc2
###################################################
print(ggs_autocorrelation(alpha0samples))


###################################################
### code chunk number 27: ggmcmc-help (eval = FALSE)
###################################################
## help(package="ggmcmc", help_type="html")


###################################################
### code chunk number 28: mcmc-DLTsamples
###################################################
DLTsamples <- mcmc(data=data1,model=DLTmodel,options=options)


###################################################
### code chunk number 29: mcmc-postDLTsamples
###################################################
data3 <-Data(x=c(25,50,50,75,100,100,225,300),
             y=c(0,0,0,0,1,1,1,1),
             doseGrid=seq(from=25,to=300,by=25))
DLTpostsamples <- mcmc(data=data3,model=DLTmodel,options=options)


###################################################
### code chunk number 30: update-DLTmodel
###################################################
newDLTmodel <- update(object=DLTmodel,data=data3)
newDLTmodel@phi1
newDLTmodel@phi2


###################################################
### code chunk number 31: plot-model-fit
###################################################
print(plot(samples, model, data))


###################################################
### code chunk number 32: empty-data
###################################################
## provide only the dose grid:
emptydata <- Data(doseGrid=data@doseGrid)
## obtain prior samples with this Data object
priorsamples <- mcmc(emptydata, model, options)
## then produce the plot
print(plot(priorsamples, model, emptydata))


###################################################
### code chunk number 33: plot-samplesdata1
###################################################
print(plot(DLTsamples,DLTmodel,data1))


###################################################
### code chunk number 34: emptydatanoDLTsamples
###################################################
print(plot(data1,DLTmodel))


###################################################
### code chunk number 35: rel-incs
###################################################
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))


###################################################
### code chunk number 36: max-dose
###################################################
nextMaxDose <- maxDose(myIncrements,
                       data=data)
nextMaxDose


###################################################
### code chunk number 37: 3fold-incs
###################################################
myIncrements1 <- IncrementsRelative(intervals=c(25),
                                   increments=c(2))


###################################################
### code chunk number 38: ncrm-spec
###################################################
myNextBest <- NextBestNCRM(target=c(0.2, 0.35),
                           overdose=c(0.35, 1),
                           maxOverdoseProb=0.25)


###################################################
### code chunk number 39: mtd-spec
###################################################
mtdNextBest <- NextBestMTD(target=0.33,
                           derive=
                               function(mtdSamples){
                                   quantile(mtdSamples, probs=0.25)
                               })


###################################################
### code chunk number 40: TD-spec
###################################################
TDNextBest <- NextBestTD(targetDuringTrial=0.35,
                         targetEndOfTrial=0.3)


###################################################
### code chunk number 41: TDsamples-spec
###################################################
TDsamplesNextBest <- NextBestTDsamples(targetDuringTrial=0.35,
                                       targetEndOfTrial=0.3,
                                       derive=function(TDsamples){
                                         quantile(TDsamples,probs=0.3)})



###################################################
### code chunk number 42: next-best-run
###################################################
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples, model=model, data=data)


###################################################
### code chunk number 43: next-best-results
###################################################
doseRecommendation$value
print(doseRecommendation$plot)


###################################################
### code chunk number 44: next-best-TD-run
###################################################
doseRecDLT <- nextBest(TDNextBest,doselimit=300,model=newDLTmodel,data=data3)


###################################################
### code chunk number 45: next-bestTD-results
###################################################
doseRecDLT$nextdose
doseRecDLT$targetDuringTrial
doseRecDLT$TDtargetDuringTrial
print(doseRecDLT$plot)


###################################################
### code chunk number 46: next-bestTDsamples-run
###################################################
doseRecDLTSamples <- nextBest(TDsamplesNextBest,doselimit=300,
                              samples=DLTpostsamples,model=newDLTmodel,
                              data=data3)


###################################################
### code chunk number 47: next-bestTDsamples-results
###################################################
print(doseRecDLTSamples$plot)


###################################################
### code chunk number 48: size-range
###################################################
mySize1 <- CohortSizeRange(intervals=c(0, 30),
                           cohortSize=c(1, 3))


###################################################
### code chunk number 49: size-dlt
###################################################
mySize2 <- CohortSizeDLT(DLTintervals=c(0, 1),
                         cohortSize=c(1, 3))


###################################################
### code chunk number 50: size-combined
###################################################
mySize <- maxSize(mySize1, mySize2)


###################################################
### code chunk number 51: size-eval
###################################################
size(mySize,
     dose=doseRecommendation$value,
     data=data)


###################################################
### code chunk number 52: size-const
###################################################
mySize <- CohortSizeConst(size=3)


###################################################
### code chunk number 53: rules-bits
###################################################
myStopping1 <- StoppingMinCohorts(nCohorts=3)
myStopping2 <- StoppingTargetProb(target=c(0.2, 0.35),
                                  prob=0.5)
myStopping3 <- StoppingMinPatients(nPatients=20)


###################################################
### code chunk number 54: rules-compose
###################################################
myStopping <- (myStopping1 & myStopping2) | myStopping3


###################################################
### code chunk number 55: rules2-bits
###################################################
myStopping4 <- StoppingTDCIRatio(targetRatio=5,
                                        targetEndOfTrial=0.3)



###################################################
### code chunk number 56: rules-try
###################################################
stopTrial(stopping=myStopping, dose=doseRecommendation$value,
          samples=samples, model=model, data=data)


###################################################
### code chunk number 57: rules2-try
###################################################
stopTrial(stopping= myStopping4, dose=doseRecDLTSamples$nextdose,
          samples=DLTpostsamples,model=newDLTmodel,data=data3)

stopTrial(stopping= myStopping4, dose=doseRecDLT$nextdose,
          model=newDLTmodel,data=data3)


###################################################
### code chunk number 58: design-setup
###################################################
design <- Design(model=model,
                 nextBest=myNextBest,
                 stopping=myStopping,
                 increments=myIncrements,
                 cohortSize=mySize,
                 data=emptydata,
                 startingDose=3)


###################################################
### code chunk number 59: TDDesign-setup
###################################################
DLTdesign <-TDDesign(model=DLTmodel,
                     nextBest=TDNextBest,
                     stopping=myStopping4,
                     increments=myIncrements1,
                     cohortSize=mySize,
                     data=data1,
                     startingDose=25)


###################################################
### code chunk number 60: TDsamplesDesign-set-up
###################################################
DLTsamplesDesign <- TDsamplesDesign(model=DLTmodel,
                                    nextBest=TDsamplesNextBest,
                                    stopping=(myStopping4|myStopping3),
                                    increments = myIncrements1,
                                    cohortSize=mySize,
                                    data=data1,
                                    startingDose=25)


###################################################
### code chunk number 61: design-examine
###################################################
set.seed(23)
examine(design)


###################################################
### code chunk number 62: true-def
###################################################
## define the true function
myTruth <- function(dose)
{
    model@prob(dose, alpha0=7, alpha1=8)
}

## plot it in the range of the dose grid
curve(myTruth(x), from=0, to=80, ylim=c(0, 1))


###################################################
### code chunk number 63: trueDLT
###################################################
## define the true function
TrueDLT <- function(dose)
{
    DLTmodel@prob(dose, phi1=-53.66584, phi2=10.50499)
}

## plot it in the range of the dose grid
curve(TrueDLT, from=25, to=300, ylim=c(0, 1))


###################################################
### code chunk number 64: run-sims
###################################################
time <- system.time(mySims <- simulate(design,
                                       args=NULL,
                                       truth=myTruth,
                                       nsim=100,
                                       seed=819,
                                       mcmcOptions=options,
                                       parallel=FALSE))[3]
time


###################################################
### code chunk number 65: sim-class
###################################################
class(mySims)


###################################################
### code chunk number 66: sim-help
###################################################
help("Simulations-class", help="html")


###################################################
### code chunk number 67: third-trial
###################################################
print(plot(mySims@data[[3]]))


###################################################
### code chunk number 68: third-dose
###################################################
mySims@doses[3]


###################################################
### code chunk number 69: third-stop
###################################################
mySims@stopReasons[[3]]


###################################################
### code chunk number 70: sim-plot
###################################################
print(plot(mySims))


###################################################
### code chunk number 71: sim-summary
###################################################
summary(mySims,
        truth=myTruth)


###################################################
### code chunk number 72: sim-sum-plot
###################################################
simSum <- summary(mySims,
                  truth=myTruth)
print(plot(simSum))


###################################################
### code chunk number 73: sim-sum-plot2
###################################################
dosePlot <- plot(simSum, type="doseSelected") +
      scale_x_continuous(breaks=10:30, limits=c(10, 30))
print(dosePlot)


###################################################
### code chunk number 74: DLTSim-run
###################################################
DLTSim <- simulate(DLTdesign,
                   args=NULL,
                   truth=TrueDLT,
                   nsim=10,
                   seed=819,
                   parallel=FALSE)


###################################################
### code chunk number 75: DLTsampSim-run
###################################################
DLTsampSim <- simulate(DLTsamplesDesign,
                       args=NULL,
                       truth=TrueDLT,
                       nsim=10,
                       seed=819, 
                       mcmcOptions=options,
                       parallel=FALSE)


###################################################
### code chunk number 76: DLTSim-dose
###################################################
DLTSim@doses[3]


###################################################
### code chunk number 77: DLTsampSim-dose
###################################################
DLTsampSim@doses[3]


###################################################
### code chunk number 78: plotDLTSim
###################################################
print(plot(DLTSim))


###################################################
### code chunk number 79: plotDLTsampSim
###################################################
print(plot(DLTsampSim))


###################################################
### code chunk number 80: DLTSim-summary
###################################################
summary(DLTSim,
        truth=TrueDLT)


###################################################
### code chunk number 81: DLTsampSim-summary
###################################################
summary(DLTsampSim,
        truth=TrueDLT)


###################################################
### code chunk number 82: DLTSim-plotsummary
###################################################
DLTsimSum <- summary(DLTSim,
             truth=TrueDLT)
print(plot(DLTsimSum))


###################################################
### code chunk number 83: DLTsampSim-plotsummary
###################################################
DLTsimsampSum <- summary(DLTsampSim,
                truth=TrueDLT)
print(plot(DLTsimsampSum))


###################################################
### code chunk number 84: explain-fut
###################################################
model@prob


###################################################
### code chunk number 85: fut-samples
###################################################
postSamples <- as.data.frame(samples@data)[(1:20)*50, ]
postSamples


###################################################
### code chunk number 86: design-future
###################################################
nowDesign <- Design(model=model,
                    nextBest=myNextBest,
                    stopping=myStopping,
                    increments=myIncrements,
                    cohortSize=mySize,
                    ## use the current data:
                    data=data,
                    ## and the recommended dose as the starting dose:
                    startingDose=doseRecommendation$value)


###################################################
### code chunk number 87: sim-future
###################################################
time <- system.time(futureSims <- simulate(
    ## supply the new design here
    nowDesign,
    ## the truth is the assumed prob function
    truth=model@prob,
    ## further arguments are the
    ## posterior samples
    args=postSamples,
    ## do exactly so many simulations as
    ## we have samples
    nsim=nrow(postSamples),
    seed=918,
    ## this remains the same:
    mcmcOptions=options,
    parallel=FALSE))[3]
time


###################################################
### code chunk number 88: sim-future-plot
###################################################
print(plot(futureSims))


###################################################
### code chunk number 89: sim-future-summary
###################################################
summary(futureSims,
        truth=myTruth)


###################################################
### code chunk number 90: three-plus-three-setup
###################################################
threeDesign <- ThreePlusThreeDesign(doseGrid=c(5, 10, 15, 25, 35, 50, 80))
class(threeDesign)


###################################################
### code chunk number 91: three-sims
###################################################
threeSims <- simulate(threeDesign,
                      nsim=1000,
                      seed=35,
                      truth=myTruth,
                      parallel=FALSE)


###################################################
### code chunk number 92: three-sims-summary
###################################################
threeSimsSum <- summary(threeSims,
                        truth=myTruth)
threeSimsSum


###################################################
### code chunk number 93: three-sims-plot
###################################################
print(plot(threeSimsSum))


###################################################
### code chunk number 94: dual-data-struct
###################################################
data <- DataDual(
    x=
        c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10,
          20, 20, 20, 40, 40, 40, 50, 50, 50),
    y=
        c(0, 0, 0, 0, 0, 0, 1, 0,
          0, 1, 1, 0, 0, 1, 0, 1, 1),
    w=
        c(0.31, 0.42, 0.59, 0.45, 0.6, 0.7, 0.55, 0.6,
          0.52, 0.54, 0.56, 0.43, 0.41, 0.39, 0.34, 0.38, 0.21),
    doseGrid=
        c(0.1, 0.5, 1.5, 3, 6,
          seq(from=10, to=80, by=2)))


###################################################
### code chunk number 95: dual-data-plot
###################################################
print(plot(data))


###################################################
### code chunk number 96: dual-rw1-model
###################################################
model <- DualEndpointRW(mu=c(0, 1),
                        Sigma=matrix(c(1, 0, 0, 1), nrow=2),
                        sigma2betaW=
                        0.01,
                        sigma2W=
                        c(a=0.1, b=0.1),
                        rho=
                        c(a=1, b=1),
                        smooth="RW1")


###################################################
### code chunk number 97: dual-options
###################################################
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=500)


###################################################
### code chunk number 98: dual-mcmc
###################################################
samples <- mcmc(data, model, options)


###################################################
### code chunk number 99: dual-conv
###################################################
data@nGrid
betaWpicks <- get(samples, "betaW", c(1, 5, 10, 25))
ggs_traceplot(betaWpicks)


###################################################
### code chunk number 100: dual-modelfit
###################################################
print(plot(samples, model, data, extrapolate=FALSE))


###################################################
### code chunk number 101: dual-variance
###################################################
ggs_histogram(get(samples, "precW"))


###################################################
### code chunk number 102: dual-nextbest
###################################################
myNextBest <- NextBestDualEndpoint(target=c(0.9,1),
                                   overdose=c(0.35, 1),
                                   maxOverdoseProb=0.25)


###################################################
### code chunk number 103: dual-nextdose-eval
###################################################
nextDose <- nextBest(myNextBest,
                     doselimit=50,
                     samples=samples,
                     model=model,
                     data=data)
nextDose$value


###################################################
### code chunk number 104: dual-nextdose-plot
###################################################
print(nextDose$plot)


###################################################
### code chunk number 105: dual-stop
###################################################
myStopping6 <- StoppingTargetBiomarker(target=c(0.9,1),
                                       prob=0.5)


###################################################
### code chunk number 106: dual-stop-try
###################################################
stopTrial(myStopping6, dose=nextDose$value,
          samples, model, data)


###################################################
### code chunk number 107: dual-stop-whole
###################################################
myStopping <- myStopping6 | StoppingMinPatients(40)


###################################################
### code chunk number 108: dual-design
###################################################
emptydata <- DataDual(doseGrid=data@doseGrid)
design <- DualDesign(model=model,
                     data=emptydata,
                     nextBest=myNextBest,
                     stopping=myStopping,
                     increments=myIncrements,
                     cohortSize=CohortSizeConst(3),
                     startingDose=6)


###################################################
### code chunk number 109: dual-scenario
###################################################
betaMod <- function (dose, e0, eMax, delta1, delta2, scal)
{
    maxDens <- (delta1^delta1) * (delta2^delta2)/((delta1 + delta2)^(delta1 + delta2))
    dose <- dose/scal
    e0 + eMax/maxDens * (dose^delta1) * (1 - dose)^delta2
}
trueBiomarker <- function(dose)
{
    betaMod(dose, e0=0.2, eMax=0.6, delta1=5, delta2=5 * 0.5 / 0.5, scal=100)
}
trueTox <- function(dose)
{
    pnorm((dose-60)/10)
}


###################################################
### code chunk number 110: dual-sc-plot
###################################################
par(mfrow=c(1, 2))
curve(trueTox(x), from=0, to=80)
curve(trueBiomarker(x), from=0, to=80)


###################################################
### code chunk number 111: dual-sims
###################################################
mySims <- simulate(design,
                   trueTox=trueTox,
                   trueBiomarker=trueBiomarker,
                   sigma2W=0.01,
                   rho=0,
                   nsim=10,
                   parallel=FALSE,
                   seed=3,
                   startingDose=6,
                   mcmcOptions =
                       McmcOptions(burnin=1000,
                                   step=1,
                                   samples=3000))


###################################################
### code chunk number 112: dual-sims-plot
###################################################
print(plot(mySims))


###################################################
### code chunk number 113: dual-sims-sum
###################################################
sumOut <- summary(mySims,
                  trueTox=trueTox,
                  trueBiomarker=trueBiomarker)
sumOut


###################################################
### code chunk number 114: dual-sim-sum-plot
###################################################
print(plot(sumOut))


###################################################
### code chunk number 115: dualdata
###################################################
data2<-DataDual(doseGrid=seq(25,300,25))

data4<-DataDual(x=c(25,50,50,75,100,100,225,300),
                y=c(0,0,0,0,1,1,1,1),
                w=c(0.31,0.42,0.59,0.45,0.6,0.7,0.6,0.52),
                doseGrid=seq(25,300,25))



###################################################
### code chunk number 116: modelEff
###################################################
Effmodel<-Effloglog(Eff=c(1.223,2.513),Effdose=c(25,300),nu=c(a=1,b=0.025),data=data2)


###################################################
### code chunk number 117: strEff
###################################################
str(Effmodel)


###################################################
### code chunk number 118: modelEffFlexi
###################################################
Effmodel2<- EffFlexi(Eff=c(1.223, 2.513),
                    Effdose=c(25,300),sigma2=c(a=0.1,b=0.1),
                    sigma2betaW=c(a=20,b=50),smooth="RW2",data=data2)


###################################################
### code chunk number 119: strEffmodel2
###################################################
str(Effmodel2)


###################################################
### code chunk number 120: mcmc-Effsamples
###################################################
Effsamples <- mcmc(data=data2,model=Effmodel,options)
Effsamples2 <- mcmc(data=data2, model=Effmodel2, options)


###################################################
### code chunk number 121: mcmc-Effpostsamples
###################################################
Effpostsamples <- mcmc(data=data2,model=Effmodel,options)
Effpostsamples2 <- mcmc(data=data2, model=Effmodel2, options)


###################################################
### code chunk number 122: update-Effmodel
###################################################
newEffmodel <- update(object=Effmodel,data=data4)
newEffmodel@theta1
newEffmodel@theta2
newEffmodel@nu


###################################################
### code chunk number 123: update-Effmodel2
###################################################
newEffmodel2 <- update(object=Effmodel2,data=data4)
newEffmodel2@RWmat


###################################################
### code chunk number 124: plot-samplesdata2loglog
###################################################
print(plot(Effpostsamples, newEffmodel,data4))


###################################################
### code chunk number 125: plot-samplesdata2Flexi
###################################################
print(plot(Effpostsamples2,newEffmodel2,data4))


###################################################
### code chunk number 126: plot-nosamplesEffmodel
###################################################
print(plot(data2,Effmodel))


###################################################
### code chunk number 127: plotDualResponseNoSamples
###################################################
plotDualResponses(DLEmodel=DLTmodel,
                  Effmodel=Effmodel,data=data2)


###################################################
### code chunk number 128: plotDualResponseSamples
###################################################
plotDualResponses(DLEmodel=DLTmodel,DLEsamples=DLTsamples,
                  Effmodel=Effmodel,Effsamples=Effsamples,data=data2)


###################################################
### code chunk number 129: plotGain
###################################################
plotGain(DLEmodel=newDLTmodel,Effmodel=newEffmodel,data=data4)


###################################################
### code chunk number 130: nextbestmaxgain
###################################################
GainNextBest <-NextBestMaxGain(DLEDuringTrialtarget=0.35,
                               DLEEndOfTrialtarget=0.3)


###################################################
### code chunk number 131: nextbest-maxgain
###################################################
doseRecGain <- nextBest(GainNextBest,
                        doselimit=max(data4@doseGrid),
                        model=newDLTmodel,
                        Effmodel=newEffmodel,
                        data=data4)


###################################################
### code chunk number 132: nextbestplot-maxgain
###################################################
doseRecGain$plot


###################################################
### code chunk number 133: nextbestmaxgainsamples
###################################################
GainsamplesNextBest <- NextBestMaxGainSamples(DLEDuringTrialtarget=0.35,
                                   DLEEndOfTrialtarget=0.3,
                                   TDderive=function(TDsamples){
                                     quantile(TDsamples,prob=0.3)},
                                   Gstarderive=function(Gstarsamples){
                                     quantile(Gstarsamples,prob=0.5)})


###################################################
### code chunk number 134: nextbest-NextBestMaxGainSamples
###################################################
doseRecGainSamples <- nextBest(GainsamplesNextBest,
                               doselimit=max(data4@doseGrid),
                               model=newDLTmodel,
                               samples=DLTpostsamples,
                               Effmodel=newEffmodel,
                               Effsamples=Effpostsamples,
                               data=data4)


###################################################
### code chunk number 135: nextbest-NextBestMaxGainSamplesplot
###################################################
doseRecGainSamples$plot


###################################################
### code chunk number 136: StoppingGstarCIratio
###################################################
myStopping7 <- StoppingGstarCIRatio(targetRatio = 5,
                                    targetEndOfTrial=0.3)

myStopping8 <- myStopping7 | StoppingMinPatients(72)



###################################################
### code chunk number 137: StopTrial-stoppingGstar
###################################################

stopTrial(stopping=myStopping7,dose=doseRecGain$nextdose,model=newDLTmodel,
          data=data4, Effmodel=newEffmodel)

stopTrial(stopping=myStopping7,
          dose=doseRecGainSamples$nextdose,
          samples=DLTpostsamples,
          model=newDLTmodel,
          data=data4,
          TDderive=function(TDsamples){
            quantile(TDsamples,prob=0.3)},
          Effmodel=newEffmodel,
          Effsamples=Effpostsamples,
          Gstarderive=function(Gstarsamples){
            quantile(Gstarsamples,prob=0.5)})


###################################################
### code chunk number 138: DualResponseDesign
###################################################
design1 <- DualResponsesDesign(nextBest=GainNextBest,
                               model=DLTmodel,
                               Effmodel=Effmodel,
                               data=data2,
                               stopping=myStopping7,
                               increments=myIncrements1,
                               cohortSize=mySize,
                               startingDose=25)

design2 <- DualResponsesSamplesDesign(nextBest=GainsamplesNextBest,
                                      model=DLTmodel,
                                      Effmodel=Effmodel,
                                      data=data2,
                                      stopping=myStopping8,
                                      increments=myIncrements1,
                                      cohortSize=mySize,
                                      startingDose=25)


###################################################
### code chunk number 139: DualResponsesDesign-Flexi
###################################################
design3 <- DualResponsesSamplesDesign(nextBest=GainsamplesNextBest,
                                      model=DLTmodel,
                                      Effmodel=Effmodel2,
                                      data=data2,
                                      stopping=myStopping8,
                                      increments=myIncrements1,
                                      cohortSize=mySize,
                                      startingDose=25)


###################################################
### code chunk number 140: trueDLTtrueEff
###################################################
myTruthDLT<- function(dose)
{ DLTmodel@prob(dose, phi1=-53.66584, phi2=10.50499)
}

myTruthEff<- function(dose)
{Effmodel@ExpEff(dose,theta1=-4.818429,theta2=3.653058)
}

myTruthGain <- function(dose)
{return(myTruthEff(dose) * (1-myTruthDLT(dose)))}


###################################################
### code chunk number 141: Truecurves
###################################################
TruthTD<-function(prob){DLTmodel@dose(prob, phi1=-53.66584, phi2=10.50499)}
GAIN<-function(xi){-(-4.8218429+3.653058*log(xi))/(1+exp(-53.66584+10.50499*xi))}
Txi<-(optim(1,GAIN,method="BFGS")$par)
maxg<-(optim(1,GAIN,method="BFGS")$value)
gstar<-exp(Txi)
td30<-TruthTD(0.3)
td35<-TruthTD(0.35)
DoseLevels<-seq(2,300,1)
plot(DoseLevels,myTruthDLT(DoseLevels), col='red',type='l',lwd=3,ylab='Values',
  ylim=c(0,max(1,max(myTruthEff(DoseLevels)))))
points(td30,0.3,col='violet',pch=15,cex=2)
points(td35,0.35,col='violet',pch=16,cex=2)
lines(DoseLevels,myTruthEff(DoseLevels),col='blue',type='l',lwd=3)
lines(DoseLevels,myTruthGain(DoseLevels),col='green3',type='l',lwd=3)
points(gstar,-maxg,col='green3',pch=17,cex=2)
legend('topright',bty='n',cex=1.2,c('p(DLT)=0.3','p(DLT)=0.35','Max gain','p(DLTs)',
'efficacy','gain'),text.col=c('violet','violet','green3','red','blue','green3'),
pch=c(15,16,17,NA,NA,NA),lty=c(NA,NA,NA,1,1,1),col=c('violet','violet','green3','red','blue','green3'))



###################################################
### code chunk number 142: TruecurveswithFlexi
###################################################
myTruthEff1<- c(-0.5478867, 0.1645417,  0.5248031,  0.7604467,  
               0.9333009  ,1.0687031,  1.1793942 , 1.2726408 , 
               1.3529598 , 1.4233411 , 1.4858613 , 1.5420182)

d1 <- data2@doseGrid
myTruthGain1 <- myTruthEff1 * (1-myTruthDLT(d1))


###################################################
### code chunk number 143: TruecurvesFlexi
###################################################
maxg1<-max(myTruthGain1)
gstar1 <- data2@doseGrid[which.max(myTruthGain1)]
DoseLevels1<-seq(1,300,1)
TruthTD<-function(prob)
  {DLTmodel@dose(prob, phi1=-53.66584, phi2=10.50499)}
td30<-TruthTD(0.3)
td35<-TruthTD(0.35)
plot(DoseLevels1,myTruthDLT(DoseLevels1), col='red',type='l',
  lwd=3,ylab='Values',ylim=c(0,max(1,max(myTruthEff1))))
points(td30,0.3,col='violet',pch=15,cex=2)
points(td35,0.35,col='violet',pch=16,cex=2)
lines(d1,myTruthEff1,col='blue',type='l',lwd=3)
lines(d1,myTruthGain1,col='green3',type='l',lwd=3)
points(gstar1,maxg1,col='green3',pch=17,cex=2)
legend('topright',bty='n',cex=1.2,c('p(DLT)=0.3','p(DLT)=0.35',
'Max gain','p(DLTs)','efficacy','gain'),text.col=c('violet','violet',
'green3','red','blue','green3'),pch=c(15,16,17,NA,NA,NA),
lty=c(NA,NA,NA,1,1,1),col=c('violet','violet','green3','red','blue','green3'))


###################################################
### code chunk number 144: simulateDualResponseDesign
###################################################
Sim1 <- simulate(object=design1,
                 args=NULL,
                 trueDLE=myTruthDLT,
                 trueEff=myTruthEff,
                 trueNu=1/0.025,
                 nsim=10,
                 seed=819,
                 parallel=FALSE)


###################################################
### code chunk number 145: simulateDualResponsesamplesDesign
###################################################
Sim2 <- simulate(object=design2,
                 args=NULL,
                 trueDLE=myTruthDLT,
                 trueEff=myTruthEff,
                 trueNu=1/0.025,
                 nsim=10,
                 seed=819,
                 mcmcOptions=options,
                 parallel=FALSE)


###################################################
### code chunk number 146: simulateDualResponseFlexiDesign
###################################################
Sim3<-simulate(object=design3,
               args=NULL,
               trueDLE=myTruthDLT,
               trueEff=myTruthEff1,
               trueSigma2=0.025,
               trueSigma2betaW=1,
               mcmcOptions=options,
               nsim=10,
               seed=819,
               parallel=FALSE)


###################################################
### code chunk number 147: plotsimulationresults1
###################################################
plot(Sim1)


###################################################
### code chunk number 148: plotsimulationresults2
###################################################
plot(Sim2)


###################################################
### code chunk number 149: plotsimulationresults3
###################################################
plot(Sim3)


###################################################
### code chunk number 150: summarysimulationresults1
###################################################
Sum1 <- summary(Sim1,
                trueDLE=myTruthDLT,
                trueEff=myTruthEff)
Sum1
print(plot(Sum1))


###################################################
### code chunk number 151: summarysimulationresults2
###################################################
Sum2 <- summary(Sim2,
                trueDLE=myTruthDLT,
                trueEff=myTruthEff)
Sum2
print(plot(Sum2))


###################################################
### code chunk number 152: summarysimulationresults3
###################################################
Sum3 <- summary(Sim3,
                trueDLE=myTruthDLT,
                trueEff=myTruthEff1)
Sum3
print(plot(Sum3))


