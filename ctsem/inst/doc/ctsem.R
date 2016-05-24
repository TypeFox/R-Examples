## ----setup, include = FALSE, cache = FALSE, echo = FALSE----------------------
library('ctsem')
library(knitr)
render_sweave()
set.seed(22)
opts_chunk$set(fig.path = 'figures/plots-', warning = FALSE, fig.align = 'center', width.cutoff = 80, fig.show = 'hold', eval = TRUE, echo = TRUE, message = FALSE, background = "white", prompt = TRUE, highlight = FALSE, comment = NA, tidy = FALSE, out.truncate = 80)
options(replace.assign = TRUE, width = 80, prompt = "R> ", scipen = 12, digits = 3)


# setwd('C:\\Users\\driver\\Dropbox\\MPIB\\CT-SEM\\manual') #set this working directory!
Sys.setenv(TEXINPUTS = getwd(),
  BIBINPUTS = getwd(),
  BSTINPUTS = getwd())

## ----eval=FALSE---------------------------------------------------------------
#  data('datastructure')
#  datastructure
#  semModel<-ctModel(n.latent=2,n.manifest=3, TRAITVAR='auto', n.TIpred=2,
#    n.TDpred=1, Tpoints=3, LAMBDA=matrix(c(1,'lambda21', 0, 0,1,0),nrow=3))
#  semFit<-ctFit(datastructure, semModel, nofit=T)
#  semFit$mxobj$A$labels
#  semFit$mxobj$S$labels
#  semFit$mxobj$M$labels
#  semFit$mxobj$F$values

## ----install, echo = T, eval = F----------------------------------------------
#  install.packages("ctsem")
#  library("ctsem")

## ----wideformat, echo = FALSE, out.truncate=100, width.cutoff=100-------------
options(width = 100)
data('datastructure')
datastructure
options(width = 80)

## ----longformat, include = TRUE, cache = FALSE, echo = FALSE, results = 'markup'----
data('longexample')
head(longexample, 7)

## ----longformatconversion, include = TRUE, cache = FALSE, echo = TRUE, results = 'markup'----
data("longexample")
wideexample <- ctLongToWide(datalong = longexample, id = "subject", 
  time = "Time", manifestNames = c("Y1", "Y2", "Y3"), 
  TDpredNames = "TD1", TIpredNames = c("TI1", "TI2"))
wide <- ctIntervalise(datawide = wideexample, Tpoints = 3, n.manifest = 3, 
  n.TDpred = 1, n.TIpred = 2, manifestNames = c("Y1", "Y2", "Y3"), 
  TDpredNames = "TD1", TIpredNames = c("TI1", "TI2") )

## ----simplemodel, include = TRUE, echo = TRUE, results = 'hide'---------------
examplemodel <- ctModel(n.latent = 2, n.manifest = 2, Tpoints = 3, 
  LAMBDA = diag(2))

## ----example1ctfit, include = TRUE, cache = TRUE, echo = TRUE, results = 'hide'----
data("ctExample1")
example1model <- ctModel(n.latent = 2, n.manifest = 2, Tpoints = 6, 
  manifestNames = c("LeisureTime", "Happiness"), 
  latentNames = c("LeisureTime", "Happiness"), LAMBDA = diag(2))
example1fit <- ctFit(datawide = ctExample1, ctmodelobj = example1model)

## ----example1ctfittable, include = TRUE, echo = TRUE--------------------------
summary(example1fit, verbose = TRUE)["discreteDRIFTstd"]

## ----expm,eval=FALSE----------------------------------------------------------
#  expm(summary(example1fit)$DRIFT * 2.5)

## ----example1testing, cache = TRUE, echo = TRUE-------------------------------
testmodel <- example1model
testmodel$DRIFT[1, 2] <- 0
testfit <- ctFit(datawide = ctExample1, ctmodelobj = testmodel)

## ----mxcompare----------------------------------------------------------------
mxCompare(example1fit$mxobj, testfit$mxobj)

## ----confidenceintervals, cache = TRUE, echo = 1------------------------------
example1cifit <- ctCI(example1fit, confidenceintervals = "DRIFT")
summary(example1cifit)$omxsummary$CI

## ----example2fit, cache = TRUE------------------------------------------------
data("ctExample1")
traitmodel <- ctModel(n.manifest = 2, n.latent = 2, Tpoints = 6, 
  LAMBDA = diag(2), manifestNames = c("LeisureTime", "Happiness"), 
  latentNames = c("LeisureTime", "Happiness"), TRAITVAR = "auto")
traitfit <- ctFit(datawide = ctExample1, ctmodelobj = traitmodel)

## ----traitparamplot, include = TRUE, cache = FALSE, echo = FALSE, results = 'hide', fig.height = 4----
par(mfrow = c(2, 2))
par(mar = c(2, 4, 2, 2))
plot(example1fit, wait = FALSE, max.time = 20, mean = FALSE, withinVariance = FALSE,randomImpulse = F,
  experimentalImpulse = F)
plot(traitfit, wait = FALSE, max.time = 20, mean = FALSE, withinVariance = FALSE,randomImpulse = F,
  experimentalImpulse = F)

## ----example1TIpred, include = TRUE, cache = TRUE, echo = TRUE, results = 'hide'----
data("ctExample1TIpred")
tipredmodel <- ctModel(n.manifest = 2, n.latent = 2, n.TIpred = 1,
  manifestNames = c("LeisureTime", "Happiness"),
  latentNames = c("LeisureTime", "Happiness"),
  TIpredNames = "NumFriends",
 Tpoints = 6, LAMBDA = diag(2), TRAITVAR = "auto")
tipredfit <- ctFit(datawide = ctExample1TIpred, ctmodelobj = tipredmodel)

summary(tipredfit, verbose = TRUE)["TIPREDEFFECT"]
summary(tipredfit, verbose = TRUE)["discreteTIPREDEFFECT"]
summary(tipredfit, verbose = TRUE)["asymTIPREDEFFECT"]
summary(tipredfit, verbose = TRUE)["addedTIPREDVAR"]

## ----example1TIpredestimates, echo = FALSE, out.width = '4cm', out.height = '4cm'----
summary(tipredfit, verbose = TRUE)["TIPREDEFFECT"]
cat("\n")
summary(tipredfit, verbose = TRUE)["discreteTIPREDEFFECT"]

## ----example1TIpredestimates2, echo = FALSE, out.width = '4cm', out.height = '4cm'----
summary(tipredfit, verbose = TRUE)["asymTIPREDEFFECT"]
cat("\n")
summary(tipredfit, verbose = TRUE)["addedTIPREDVAR"]

## ----tdpreddemo, echo = FALSE, eval = TRUE, fig.height = 3--------------------
Tpoints = 20
testm <- ctModel(Tpoints = Tpoints, n.latent = 1, n.TDpred = 1, n.manifest = 1, LAMBDA = diag(1), DRIFT = diag(-.3, 1),
  DIFFUSION = diag(.1, 1), TDPREDEFFECT = diag(1, 1), TDPREDVAR = diag(0, Tpoints-1), CINT = diag(4, 1), T0MEANS = matrix(0, ncol = 1, nrow = 1),
  T0VAR = diag(100, 1), TDPREDMEANS = matrix(c(0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ncol = 1, nrow = (Tpoints-1)))
testd <- ctGenerate(testm, n.subjects = 100, burnin = 300)

testm <- ctModel(Tpoints = Tpoints, n.latent = 1, n.TDpred = 1, n.manifest = 1, LAMBDA = diag(1), DRIFT = diag(-.3, 1),
  DIFFUSION = diag(.0, 1), TDPREDEFFECT = diag(1, 1), TDPREDVAR = diag(0, Tpoints-1), CINT = diag(4, 1), T0MEANS = matrix(0, ncol = 1, nrow = 1),
  T0VAR = diag(1, 1), TDPREDMEANS = matrix(c(0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ncol = 1, nrow = (Tpoints-1)))
testdpure <- ctGenerate(testm, n.subjects = 1, burnin = 300)

par(mfrow = c(1, 2), cex = .7)
plot(0:(Tpoints-1), testdpure[, 1:Tpoints], type = 'l', ylim = c(min(testd[, 1:Tpoints]), max(testd[, 1:Tpoints])),
  ylab = 'Dependent variable', xlab = 'Time', lwd = 3, main = 'Impulse predictor')
for(i in 1:5){
  points(0:(Tpoints-1), testd[i, 1:Tpoints], col = 1+i, type = 'b')
}

Tpoints = 20
testm <- ctModel(Tpoints = Tpoints, n.latent = 1, n.TDpred = 1, n.manifest = 1, LAMBDA = diag(1), DRIFT = diag(-.3, 1),
  DIFFUSION = diag(.15, 1), TDPREDEFFECT = diag(1.6, 1), TDPREDVAR = diag(0, Tpoints-1), CINT = diag(4, 1), T0MEANS = matrix(0, ncol = 1, nrow = 1),
  T0VAR = diag(100, 1), TDPREDMEANS = matrix(c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), ncol = 1, nrow = (19)))
testd <- ctGenerate(testm, n.subjects = 100, burnin = 300)

testm <- ctModel(Tpoints = Tpoints, n.latent = 1, n.TDpred = 1, n.manifest = 1, LAMBDA = diag(1), DRIFT = diag(-.3, 1),
  DIFFUSION = diag(.0, 1), TDPREDEFFECT = diag(1.6, 1), TDPREDVAR = diag(0, Tpoints-1), CINT = diag(4, 1), T0MEANS = matrix(0, ncol = 1, nrow = 1),
  T0VAR = diag(1, 1), TDPREDMEANS = matrix(c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), ncol = 1, nrow = (19)))
testdpure <- ctGenerate(testm, n.subjects = 1, burnin = 300)

plot(0:(Tpoints-1), testdpure[, 1:Tpoints], type = 'l', ylim = c(min(testd[, 1:Tpoints]), max(testd[, 1:Tpoints])), ylab = 'Dependent variable', xlab = 'Time', lwd = 3, main = 'Level predictor')
for(i in 1:5){
  points(0:(Tpoints-1), testd[i, 1:Tpoints], col = 1+i, type = 'b')
}

## ----example2TDpred, include = TRUE, cache = TRUE, echo = TRUE, results = 'hide', fig.height = 4, fig.width = 4, fig.align = 'center'----
data("ctExample2")
tdpredmodel <- ctModel(n.manifest = 2, n.latent = 2, n.TDpred = 1, 
  Tpoints = 8, manifestNames = c("LeisureTime", "Happiness"), 
  TDpredNames = "MoneyInt", latentNames = c("LeisureTime", "Happiness"),
  T0TDPREDCOV = matrix(0, nrow = 2, ncol=7),
  TRAITTDPREDCOV = matrix(0, nrow = 2, ncol=7), 
  LAMBDA = diag(2), TRAITVAR = "auto")
tdpredfit <- ctFit(datawide = ctExample2, ctmodelobj = tdpredmodel)

summary(tdpredfit, verbose = TRUE)["TDPREDEFFECT"]
summary(tdpredfit, verbose = TRUE)["discreteTDPREDEFFECT"]

## ----example2TDpredestimates, echo = FALSE, out.width = '3cm'-----------------
summary(tdpredfit, verbose = TRUE)["TDPREDEFFECT"]

## ----example2TDpredestimates3, echo = FALSE, out.width = '3cm'----------------
summary(tdpredfit, verbose = TRUE)["discreteTDPREDEFFECT"]

## ----example2TDpredlevel, include = TRUE, eval = TRUE, echo = TRUE, cache = TRUE----
data("ctExample2")
tdpredlevelmodel <- ctModel(n.manifest = 2, n.latent = 3, n.TDpred = 1, 
  Tpoints = 8, manifestNames = c("LeisureTime", "Happiness"), 
  TDpredNames = "MoneyInt", 
  latentNames = c("LeisureTime", "Happiness", "MoneyIntLatent"),
  T0TDPREDCOV = matrix(0, nrow = 3, ncol = 7),
  TRAITTDPREDCOV = matrix(0, nrow = 3, ncol = 7), 
  LAMBDA = matrix(c(1,0, 0,1, 0,0), ncol = 3), TRAITVAR = "auto")

tdpredlevelmodel$TRAITVAR[3, ] <- 0
tdpredlevelmodel$TRAITVAR[, 3] <- 0
tdpredlevelmodel$DIFFUSION[, 3] <- 0
tdpredlevelmodel$DIFFUSION[3, ] <- 0
tdpredlevelmodel$T0VAR[3, ] <- 0
tdpredlevelmodel$T0VAR[, 3] <- 0
tdpredlevelmodel$CINT[3] <- 0
tdpredlevelmodel$T0MEANS[3] <- 0
tdpredlevelmodel$TDPREDEFFECT[3, ] <- 1
tdpredlevelmodel$DRIFT[3, ] <- 0

tdpredlevelfit <- ctFit(datawide = ctExample2, 
  ctmodelobj = tdpredlevelmodel)

summary(tdpredlevelfit, verbose = TRUE)["DRIFT"]
summary(tdpredlevelfit, timeInterval = 20, 
  verbose = TRUE)["discreteTDPREDEFFECT"]

## ----timeseries, cache = TRUE, echo = TRUE------------------------------------
data("ctExample3")
model <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 100, 
  LAMBDA = matrix(c(1, "lambda2", "lambda3"), nrow = 3, ncol = 1), 
  MANIFESTMEANS = matrix(c(0, "manifestmean2", "manifestmean3"), nrow = 3, 
    ncol = 1))
fit <- ctFit(data = ctExample3, ctmodelobj = model, objective = "Kalman",
  stationary = c("T0VAR"))

## ----multigroup, cache = TRUE, echo = TRUE------------------------------------
data("ctExample4")

basemodel <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 20,
  LAMBDA = matrix(c(1, "lambda2", "lambda3"), nrow = 3, ncol = 1),
  TRAITVAR="auto", MANIFESTMEANS = matrix(c(0, "manifestmean2", 
    "manifestmean3"), nrow = 3, ncol = 1))

freemodel <- basemodel
freemodel$LAMBDA[3, 1] <- "groupfree"
groups <- paste0("g", rep(1:2, each = 10))

multif <- ctMultigroupFit(datawide = ctExample4, groupings = groups,
  ctmodelobj = basemodel, freemodel = freemodel)

## ----multigroupOutput, echo = FALSE-------------------------------------------
multif$mxobj$output$estimate[grep("lambda3", names(multif$mxobj$output$estimate))]

## ----dynamicresidualmodel, echo = TRUE, results='hide',fig.keep='none'--------
genm <- ctModel(Tpoints = 200, n.latent = 2, n.manifest = 1, 
  LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
  DIFFUSION = matrix(c(0, 0, 0, 1), 2),
  MANIFESTVAR = diag(.6,1),
  DRIFT = matrix(c(0, -.1, 1, -.2), nrow = 2),   
  CINT = matrix(c(1, 0), nrow = 2))

data <- ctGenerate(genm, n.subjects = 1, burnin = 200, dT = 1)

ctIndplot(data, n.subjects = 1 , n.manifest = 1, Tpoints = 200)

model <- ctModel(Tpoints = 200, n.latent = 2, n.manifest = 1, 
  LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
  DIFFUSION = matrix(c(0, 0, 0, "diffusion"), 2),
  DRIFT = matrix(c(0, "regulation", 1, "diffusionAR"), nrow = 2),   
  CINT = matrix(c("processCINT", 0), nrow = 2))

fit <- ctFit(data, model, stationary = c("T0VAR"), carefulFit = FALSE)

## ----osscilating, cache = TRUE, echo = TRUE, eval = TRUE, include = TRUE------
data("Oscillating")

inits <- c(-38, -.5, 1, 1, .1, 1, 0, .9)
names(inits) <- c("crosseffect","autoeffect", "diffusion",
  "T0var11", "T0var21", "T0var22","m1", "m2")

oscillatingm <- ctModel(n.latent = 2, n.manifest = 1, Tpoints = 11, 
  MANIFESTVAR = matrix(c(0), nrow = 1, ncol = 1), 
  LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
  T0MEANS = matrix(c('m1', 'm2'), nrow = 2, ncol = 1), 
  T0VAR = matrix(c("T0var11", "T0var21", 0, "T0var22"), nrow = 2, ncol = 2),
  DRIFT = matrix(c(0, "crosseffect", 1, "autoeffect"), nrow = 2, ncol = 2), 
  CINT = matrix(0, ncol = 1, nrow = 2),
  DIFFUSION = matrix(c(0, 0, 0, "diffusion"), nrow = 2, ncol = 2),
  startValues = inits)

oscillatingf <- ctFit(Oscillating, oscillatingm, carefulFit = FALSE)

## ----inits, echo = TRUE-------------------------------------------------------
omxInits <- omxGetParameters(example1fit$mxobj)

fitWithInits <- ctFit(data = ctExample1, ctmodelobj = example1model, 
  omxStartValues = omxInits)

## ----fulldiscretecomparison,cache=TRUE----------------------------------------
data("ctExample1")
traitmodel <- ctModel(n.manifest = 2, n.latent = 2, Tpoints = 6, 
  LAMBDA = diag(2), manifestNames = c("LeisureTime", "Happiness"), 
  latentNames = c("LeisureTime", "Happiness"), TRAITVAR = "auto")
traitfit <- ctFit(datawide = ctExample1, ctmodelobj = traitmodel)

discrete <- ctExample1
discrete[ , paste0('dT', 1:5)] <- mean(discrete[ , paste0('dT', 1:5)])
discretefit <- ctFit(discrete, traitmodel)

discretefit <- ctCI(discretefit, confidenceintervals = 'DRIFT')
summary(discretefit)$omxsummary$CI

## ----traitfitconfidenceintervals,cache=TRUE-----------------------------------
example1traitfit <- ctCI(traitfit, confidenceintervals = 'DRIFT')
summary(example1traitfit)$omxsummary$CI

## ----packagecomparison, eval=FALSE, cache=TRUE, fig.keep='none'---------------
#    if (!requireNamespace("PSM", quietly = TRUE)) {
#      stop("PSM package needed for this function to work. Please install it.",
#        call. = FALSE)
#    }
#    if (!requireNamespace("cts", quietly = TRUE)) {
#      stop("cts package needed for this function to work. Please install it.",
#        call. = FALSE)
#    }
#    if (!requireNamespace("yuima", quietly = TRUE)) {
#      stop("yuima package needed for this function to work. Please install it.",
#        call. = FALSE)
#    }
#  output <- matrix(NA, 10, 7)
#  colnames(output) <- c('True', 'ctsem', 'cts', 'PSM', 'yuima', 'arima', 'OpenMx')
#  for(i in 1:nrow(output)){
#  generatingModel <- ctModel(n.latent = 1, n.manifest = 1, Tpoints = 500,
#    LAMBDA = diag(1), DRIFT = matrix(-.3, nrow = 1),
#    CINT = matrix(3, 1, 1),
#    MANIFESTVAR = diag(0, 1),
#    DIFFUSION = t(chol(diag(5, 1))))
#  
#  output[i, 1] <- generatingModel$DRIFT #true value
#  
#  ctsemData <- ctGenerate(generatingModel, n.subjects=1, burnin=300)
#  longData <- ctWideToLong(ctsemData, Tpoints=500, n.manifest=1)
#  longData <- ctDeintervalise(longData)
#  
#  ### ctsem package
#  ctsemModel <- ctModel(n.latent=1, n.manifest = 1, Tpoints = 500,
#    MANIFESTVAR = diag(0, 1),
#    LAMBDA = diag(1))
#  ctsemFit <- ctFit(ctsemData, ctsemModel, stationary = c('T0VAR'))
#  output[i,2] <- mxEval(DRIFT, ctsemFit$mxobj)
#  
#  ### CTS package
#  ctsData <- longData[, c('AbsTime', 'Y1')]
#  library(cts)
#  ctsFit <- car(ctsData, order = 1, scale = 1)
#  output[i, 3] <- -1 * (1 + ctsFit$phi) / (1 - ctsFit$phi)
#  
#  ### PSM package
#  psmFit <- ctPSMfit(ctsemData, omxStartValues =
#      omxGetParameters(ctsemFit$mxobj), ctsemModel)
#  output[i, 4] <- -exp(psmFit$PSMfit$opt$par[2])
#  
#  ### yuima package (not plotted - potential issues due to dT=1)
#  library(yuima)
#  mod <- setModel(drift="drift * x + cint", diffusion = "diffusion")
#  ou <- setYuima(model = mod, data = setData(longData[,'Y1'], delta = 1))
#  mlout <- qmle(ou,start = list(drift = -.3, diffusion = 1, cint = 1))
#  output[i, 5] <- mlout@coef[2]
#  
#  ### arima (from stats package, discrete time analysis only)
#  arfit <- arima(longData[, 'Y1'], order = c(1, 0, 0))
#  log(arfit$coef[1]) #transform ar1 parameter to drift parameter
#  output[i, 6]<-log(arfit$coef[1])
#  
#  ### OpenMx state space continuous time function (specified via ctsem here)
#  ctsemModel <- ctModel(n.latent=1, n.manifest = 1, Tpoints = 500,
#    MANIFESTVAR = diag(0, 1), T0VAR=diag(1), LAMBDA = diag(1))
#  mxFit <- ctFit(ctsemData, ctsemModel, objective='Kalmanmx',
#    carefulFit = FALSE)
#  output[i,7] <- mxEval(DRIFT, ctsemFit$mxobj)
#  
#  } #end for loop
#  
#  ### plot output
#  plot(density(output[1:i, 2]), xlim = c(-.4, -.1), lty = 2, lwd = 1,
#    ylab='Density',
#    main = 'Density of estimates of drift parameter (true value -0.3)')
#  points(density(output[1:i, 3]), col='red', type='l', lty=3, lwd=1)
#  points(density(output[1:i, 4]), col='green', type='l', lty=2, lwd=1)
#  points(density(output[1:i, 6]), col='blue', type='l', lwd=1, lty=3)
#  points(density(output[1:i, 7]), col='yellow', type='l', lwd=1, lty=5)
#  legend('topleft', bty='n',
#    text.col = c('black','red','green','blue', 'yellow'),
#    legend = c('ctsem', 'cts', 'PSM', 'arima', 'OpenMx'))

