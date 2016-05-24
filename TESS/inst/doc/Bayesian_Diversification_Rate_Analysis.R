## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(concordance=TRUE)
knit_theme$set("default")

if ( file.exists("results") == FALSE ) dir.create("results")


## ----echo=TRUE, message=FALSE--------------------------------------------
library(TESS)     # load the package
data(conifers)    # load the conifers dataset

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  myTree <- read.nexus("data/myTree.nex")

## ----echo=TRUE-----------------------------------------------------------
times <- as.numeric( branching.times(conifers) )

## ----echo=TRUE, label=plotConifers, include=TRUE, fig.cap="Conifer phylogeny from \\cite{Leslie2012} without taxon labels.", fig.pos='!ht'----
plot(conifers,show.tip.label=FALSE,no.margin=TRUE)

## ----echo=TRUE, label=lttConifers, include=TRUE, fig.cap='Lineage-through-time plot of the conifer phylogeny.',fig.pos='!ht'----
ltt.plot(conifers,log="y")

## ----echo = TRUE, eval=FALSE---------------------------------------------
#  speciation <- 1.0
#  extinction <- 0.0
#  tmrca <- 3.0

## ----echo = TRUE, eval=FALSE---------------------------------------------
#  trees <- tess.sim.age(n = 50,
#                        age = tmrca,
#                        lambda = speciation,
#                        mu = extinction,
#                        MRCA = TRUE)

## ----echo = TRUE, eval=FALSE---------------------------------------------
#  mltt.plot(trees,
#            log = "y",
#            dcol = FALSE,
#            legend = FALSE,
#            backward = FALSE)

## ----echo = TRUE, include = TRUE, eval=FALSE-----------------------------
#  expected <- function(t)
#    tess.nTaxa.expected(begin = 0,
#                        t = t,
#                        end = tmrca,
#                        lambda = speciation,
#                        mu = extinction,
#                        MRCA = TRUE,
#                        reconstructed = TRUE)
#  
#  curve(expected,add=FALSE,col="red",lty=2,lwd=5)
#  legend("topleft",col="red",lty=2,"Expected Diversity")

## ----echo=TRUE-----------------------------------------------------------
speciation <- 5.0
extinction <- 4.0
tmrca <- 3.0

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  trees <- tess.sim.age(n = 50,
#                        age = tmrca,
#                        lambda = speciation,
#                        mu = extinction,
#                        MRCA = TRUE)

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  mltt.plot(trees,
#            log = "y",
#            dcol = FALSE,
#            legend = FALSE,
#            backward = FALSE)

## ----echo=TRUE, eval = FALSE---------------------------------------------
#  expected <- function(t)
#    tess.nTaxa.expected(begin = 0,
#                        t = t,
#                        end = tmrca,
#                        lambda = speciation,
#                        mu = extinction,
#                        MRCA = TRUE,
#                        reconstructed = TRUE)
#  
#  curve(expected,add=TRUE,col="red",lty=2,lwd=5)
#  legend("topleft",col="red",lty=2,"Expected Diversity")

## ----echo=TRUE-----------------------------------------------------------
speciation <- function(t) 0.5 + 2 * exp(-1.0*t)
extinction <- 0.0
tmrca <- 3.0

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  trees <- tess.sim.age(n = 50,
#                        age = tmrca,
#                        lambda = speciation,
#                        mu = extinction,
#                        MRCA = TRUE)

## ----echo=TRUE, eval = FALSE---------------------------------------------
#  mltt.plot(trees,
#            log = "y",
#            dcol = FALSE,
#            legend = FALSE,
#            backward = FALSE)

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  expected <- function(t)
#    tess.nTaxa.expected(begin = 0,
#                        t = t,
#                        end = tmrca,
#                        lambda = speciation,
#                        mu = extinction,
#                        MRCA = TRUE,
#                        reconstructed = TRUE)
#  
#  curve(expected,add=TRUE,col="red",lty=2,lwd=5)
#  legend("topleft",col="red",lty=2,"Expected Diversity")

## ----echo=FALSE, label=lttPlots, include=TRUE, fig.height=3, out.width="\\linewidth", fig.align="center", fig.cap="Lineage-through-time curves for pure-birth trees (panel A), birth-death trees (panel B), and pure-birth trees with exponentially decreasing speciation rate (panel C)."----

# Pure-birth
pureBirthSpeciation <- 1.0
pureBirthExtinction <- 0.0
tmrca <- 3.0

pureBirthTrees <- tess.sim.age(n = 50,
                                age = tmrca,
                                lambda = pureBirthSpeciation,
                                mu = pureBirthExtinction,
                                MRCA = TRUE)


# Birth-death
birthDeathSpeciation <- 5.0
birthDeathExtinction <- 4.0
tmrca <- 3.0

birthDeathTrees <- tess.sim.age(n = 50,
                                age = tmrca,
                                lambda = birthDeathSpeciation,
                                mu = birthDeathExtinction,
                                MRCA = TRUE)


# Decreasing rate pure-birth trees
rateDecreaseSpeciation <- function(t) 0.5 + 2 * exp(-1.0*t)
rateDecreaseExtinction <- 0.0
tmrca <- 3.0

rateDecreaseTrees <- tess.sim.age(n = 50,
                                  age = tmrca,
                                  lambda = rateDecreaseSpeciation,
                                  mu = rateDecreaseExtinction,
                                  MRCA = TRUE)

# Expected number of taxa function
expected <- function(t,speciation,extinction)
  tess.nTaxa.expected(begin = 0,
                            t = t,
                            end = tmrca,
                            lambda = speciation,
                            mu = extinction,
                            MRCA = TRUE,
                            reconstructed = TRUE)

par(mfrow=c(1,3),mar=c(5,4,3,0.1),las=1)

# Plot the trees
mltt.plot(pureBirthTrees,log = "y",dcol = FALSE, legend = FALSE,backward = FALSE)
mtext("A", line = 1)
curve(expected(t=x,speciation=pureBirthSpeciation,extinction=pureBirthExtinction),add=TRUE,col="red",lty=2,lwd=2)
legend("topleft",col="red",lty=2,"Expected Diversity",lwd=2,bty='n')

mltt.plot(birthDeathTrees,log = "y",dcol = FALSE,legend = FALSE,backward = FALSE)
mtext("B", line = 1)
curve(expected(t=x,speciation=birthDeathSpeciation,extinction=birthDeathExtinction),add=TRUE,col="red",lty=2,lwd=2)
legend("topleft",col="red",lty=2,"Expected Diversity",lwd=2,bty='n')

mltt.plot(rateDecreaseTrees,log = "y",dcol = FALSE,legend = FALSE,backward = FALSE)
mtext("C", line = 1)
curve(expected(t=x,speciation=rateDecreaseSpeciation,extinction=rateDecreaseExtinction),add=TRUE,col="red",lty=2,lwd=2)
legend("topleft",col="red",lty=2,"Expected Diversity",lwd=2,bty='n')


## ----echo = TRUE---------------------------------------------------------
prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
priorsConstBD <- c("diversification"=prior_delta,
                   "turnover"=prior_tau)

## ----echo=TRUE-----------------------------------------------------------
likelihoodConstBD <- function(params) {

  speciation <- params[1] + params[2]
  extinction <- params[2]

  lnl <- tess.likelihood(times,
                         lambda = speciation,
                         mu = extinction,
                         samplingProbability = 1.0,
                         log = TRUE)

  return (lnl)

}

## ----echo = TRUE, eval=FALSE---------------------------------------------
#  set.seed(12345)   # remove this line to obtain a random seed
#  samplesConstBD <- tess.mcmc(likelihoodFunction = likelihoodConstBD,
#                              priors = priorsConstBD,
#                              parameters = runif(2,0,1),
#                              logTransforms = c(TRUE,TRUE),
#                              delta = c(1,1),
#                              iterations = 10000,
#                              burnin = 1000,
#                              thinning = 10,
#                              adaptive = TRUE,
#                              verbose = TRUE)

## ----echo = FALSE, eval=TRUE---------------------------------------------
if ( file.exists("results/samplesConstBD.rds") == FALSE ) {
set.seed(12345)   # remove this line to obtain a random seed
samplesConstBD <- tess.mcmc(likelihoodFunction = likelihoodConstBD,
                            priors = priorsConstBD,
                            parameters = runif(2,0,1),
                            logTransforms = c(TRUE,TRUE),
                            delta = c(1,1),
                            iterations = 10000,
                            burnin = 1000,
                            thinning = 10,
                            adaptive = TRUE,
                            verbose = TRUE)
saveRDS(samplesConstBD, "results/samplesConstBD.rds")
} else {
samplesConstBD <- readRDS("results/samplesConstBD.rds")
}

## ----echo = TRUE---------------------------------------------------------
summary(samplesConstBD)

## ----label=mcmcConstBD, include=TRUE, fig.cap="Trace plots (left) and marginal posterior probability densities (right) for the diversification rate (top) and turnover rate (bottom) from the MCMC simulation under the constant-rate birth-death process.", fig.pos='!ht'----
plot(samplesConstBD)

## ----echo = TRUE---------------------------------------------------------
prior_delta <- function(x) { dexp(x,rate=0.1,log=TRUE) }
prior_lambda <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_alpha <- function(x) { dexp(x,rate=0.1,log=TRUE) }
priorsDecrBD <- c("turnover"=prior_delta,
                  "initial speciation"=prior_lambda,
                  "speciation decay"=prior_alpha)

## ----echo=TRUE-----------------------------------------------------------
likelihoodDecrBD <- function(params) {

  speciation <- function(t) params[1] + params[2] * exp(-params[3]*t)
  extinction <- function(t) params[1]

  lnl <- tess.likelihood(times,
                         lambda = speciation,
                         mu = extinction,
                         samplingProbability = 1.0,
                         log = TRUE)

  return (lnl)

}

## ----echo = TRUE, eval=FALSE---------------------------------------------
#  set.seed(12345)
#  samplesDecrBD <- tess.mcmc(likelihoodFunction = likelihoodDecrBD,
#                             priors = priorsDecrBD,
#                             parameters = runif(3,0,1),
#                             logTransforms = c(TRUE,TRUE,TRUE),
#                             delta = c(1,1,1),
#                             iterations = 10000,
#                             burnin = 1000,
#                             thinning = 10,
#                             adaptive = TRUE,
#                             verbose = TRUE)

## ----echo = FALSE, eval=TRUE---------------------------------------------
if ( file.exists("results/samplesDecrBD.rds") == FALSE ) {
set.seed(12345)
samplesDecrBD <- tess.mcmc(likelihoodFunction = likelihoodDecrBD,
                           priors = priorsDecrBD,
                           parameters = runif(3,0,1),
                           logTransforms = c(TRUE,TRUE,TRUE),
                           delta = c(1,1,1),
                           iterations = 10000,
                           burnin = 1000,
                           thinning = 10,
                           adaptive = TRUE,
                           verbose = TRUE)
saveRDS(samplesDecrBD, "results/samplesDecrBD.rds")
} else {
samplesDecrBD <- readRDS("results/samplesDecrBD.rds")
}

## ----echo = TRUE---------------------------------------------------------
summary(samplesDecrBD)

## ----label=mcmcPlotDecrBD, echo = TRUE, include=TRUE, fig.cap="Trace plots and estimated posterior distribution of the parameter under the decreasing speciation rate birth-death model.", eval=TRUE----
plot(samplesDecrBD)

## ----echo = TRUE---------------------------------------------------------
rateChangeTime <- max( times ) / 2

## ----echo = TRUE---------------------------------------------------------
prior_delta_before <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau_before <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_delta_after <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau_after <- function(x) { dexp(x,rate=10.0,log=TRUE) }
priorsEpisodicBD <- c("diversification before"=prior_delta_before,
                      "turnover before"=prior_tau_before,
                      "diversification after"=prior_delta_after,
                      "turnover after"=prior_tau_after)

## ----echo=TRUE-----------------------------------------------------------
likelihoodEpisodicBD <- function(params) {

  speciation <- c(params[1]+params[2],params[3]+params[4])
  extinction <- c(params[2],params[4])

  lnl <- tess.likelihood.rateshift(times,
                                   lambda = speciation,
                                   mu = extinction,
                                   rateChangeTimesLambda = rateChangeTime,
                                   rateChangeTimesMu = rateChangeTime,
                                   samplingProbability = 1.0,
                                   log = TRUE)

  return (lnl)

}

## ----echo = TRUE, eval=FALSE---------------------------------------------
#  set.seed(12345)
#  samplesEpisodicBD <- tess.mcmc(likelihoodFunction = likelihoodEpisodicBD,
#                                 priors = priorsEpisodicBD,
#                                 parameters = runif(4,0,1),
#                                 logTransforms = c(TRUE,TRUE,TRUE,TRUE),
#                                 delta = c(1,1,1,1),
#                                 iterations = 10000,
#                                 burnin = 1000,
#                                 thinning = 10,
#                                 adaptive = TRUE,
#                                 verbose = TRUE)

## ----echo = FALSE, eval=TRUE---------------------------------------------
if ( file.exists("results/samplesEpisodicBD.rds") == FALSE ) {
set.seed(12345)
samplesEpisodicBD <- tess.mcmc(likelihoodFunction = likelihoodEpisodicBD,
                               priors = priorsEpisodicBD,
                               parameters = runif(4,0,1),
                               logTransforms = c(TRUE,TRUE,TRUE,TRUE),
                               delta = c(1,1,1,1),
                               iterations = 10000,
                               burnin = 1000,
                               thinning = 10,
                               adaptive = TRUE,
                               verbose = TRUE)
saveRDS(samplesEpisodicBD, "results/samplesEpisodicBD.rds")
} else {
samplesEpisodicBD <- readRDS("results/samplesEpisodicBD.rds")
}

## ----echo = TRUE---------------------------------------------------------
summary(samplesEpisodicBD)

## ----label=mcmcPlotEpisodicBD, echo = TRUE, include=TRUE, fig.cap="Trace plots and estimated posterior distributions of the parameters under a birth-death-shift model.", eval=TRUE----
plot(samplesEpisodicBD)

## ----echo = TRUE---------------------------------------------------------
survivalProbability <- 0.1

## ----echo = TRUE---------------------------------------------------------
prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_time <- function(x) { dunif(x,min=max(times)/2,max=max(times),log=TRUE)}
priorsMassExtinctionBD <- c("diversification"=prior_delta,
                            "turnover"=prior_tau,
                            "mass-extinction time"=prior_time)

## ----echo=TRUE-----------------------------------------------------------
likelihoodMassExtinctionBD <- function(params) {

  speciation <- params[1]+params[2]
  extinction <- params[2]
  time <- params[3]

  lnl <- tess.likelihood(times,
                         lambda = speciation,
                         mu = extinction,
                         massExtinctionTimes = time,
                         massExtinctionSurvivalProbabilities =
                                      survivalProbability,
                         samplingProbability = 1.0,
                         log = TRUE)

  return (lnl)

}

## ----echo = TRUE, eval=FALSE---------------------------------------------
#  set.seed(12345)
#  samplesMassExtinctionBD <- tess.mcmc(likelihoodFunction =
#                                         likelihoodMassExtinctionBD,
#                                         priors = priorsMassExtinctionBD,
#                                         parameters = c(runif(2,0,1),max(times)*3/4),
#                                         logTransforms = c(TRUE,TRUE,FALSE),
#                                         delta = c(1,1,1),
#                                         iterations = 10000,
#                                         burnin = 1000,
#                                         thinning = 10,
#                                         adaptive = TRUE,
#                                         verbose = TRUE)

## ----echo = FALSE, eval=TRUE---------------------------------------------
if ( file.exists("results/samplesMassExtinctionBD.rds") == FALSE ) {
set.seed(12345)
samplesMassExtinctionBD <- tess.mcmc(likelihoodFunction =
                                       likelihoodMassExtinctionBD,
                                       priors = priorsMassExtinctionBD,
                                       parameters = c(runif(2,0,1),max(times)*3/4),
                                       logTransforms = c(TRUE,TRUE,FALSE),
                                       delta = c(1,1,1),
                                       iterations = 10000,
                                       burnin = 1000,
                                       thinning = 10,
                                       adaptive = TRUE,
                                       verbose = TRUE)
saveRDS(samplesMassExtinctionBD, "results/samplesMassExtinctionBD.rds")
} else {
samplesMassExtinctionBD <- readRDS("results/samplesMassExtinctionBD.rds")
}

## ----echo = TRUE---------------------------------------------------------
summary(samplesMassExtinctionBD)

## ----label=mcmcPlotMassExtinctionBD, echo = TRUE, include=TRUE, eval=TRUE----
plot(samplesMassExtinctionBD)

## ----echo=TRUE, label=lttPlotsSampling, include=TRUE, fig.height=3, out.width="\\linewidth", fig.align="center", fig.cap="Lineage-through-time plots for completely sampled trees (panel A), incomplete trees with uniform sampling (panel B), and incomplete trees with diversified sampling (panel C)."----

# Birth-death
birthDeathSpeciationSampling <- 2.0
birthDeathExtinctionSampling <- 1.0

birthDeathTreesComplete <- tess.sim.age(n = 50,
                                age = 3.0,
                                lambda = birthDeathSpeciationSampling,
                                mu = birthDeathExtinctionSampling,
                                MRCA = TRUE)

birthDeathTreesUniform <- tess.sim.age(n = 50,
                                age = 4.0,
                                lambda = birthDeathSpeciationSampling,
                                mu = birthDeathExtinctionSampling,
                                samplingProbability = 0.25,
                                samplingStrategy = "uniform",
                                MRCA = TRUE)

birthDeathTreesDiversified <- tess.sim.age(n = 50,
                                age = 4.0,
                                lambda = birthDeathSpeciationSampling,
                                mu = birthDeathExtinctionSampling,
                                samplingProbability = 0.25,
                                samplingStrategy = "diversified",
                                MRCA = TRUE)



par(mfrow=c(1,3),mar=c(5,4,3,0.1),las=1)

# Plot the trees
mltt.plot(birthDeathTreesComplete,log = "y",dcol = FALSE,
          legend = FALSE,backward = FALSE)
mtext("A", line = 1)

mltt.plot(birthDeathTreesUniform,log = "y",dcol = FALSE,
          legend = FALSE,backward = FALSE)
mtext("B", line = 1)

mltt.plot(birthDeathTreesDiversified,log = "y",dcol = FALSE,
          legend = FALSE,backward = FALSE)
mtext("C", line = 1)


## ----echo=TRUE, include=TRUE---------------------------------------------
# simulate a tree under diversified taxon sampling
tree.diversified <- tess.sim.age(n = 1,
                                age = 4.0,
                                lambda = 2.0,
                                mu = 1.0,
                                samplingProbability = 0.25,
                                samplingStrategy = "diversified",
                                MRCA = TRUE)[[1]]

# extract the branching times from the tree
times.diversified <- as.numeric( branching.times(tree.diversified) )

## ----label=treeSamplingDiversified, include=TRUE, fig.height=4, fig.width=4, fig.align='center', fig.cap="The simulated tree under \\emph{diversified} sampling with sampling fraction $\\rho=0.25$.", fig.pos='!ht'----
plot(tree.diversified,show.tip.label=FALSE,no.margin=TRUE)

## ----echo = TRUE---------------------------------------------------------
# There are 630 known conifer species
samplingFraction <- (conifers$Nnode + 1) / 630

## ----echo = TRUE---------------------------------------------------------
prior_delta <- function(x) { dexp(x,rate=1.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=1.0,log=TRUE) }
priorsSampling <- c("diversifiation"=prior_delta,
                    "turnover"=prior_tau)

## ----echo=TRUE-----------------------------------------------------------
likelihoodUniform <- function(params) {

  speciation <- params[1] + params[2]
  extinction <- params[2]

  lnl <- tess.likelihood(times.diversified,
                         lambda = speciation,
                         mu = extinction,
                         samplingProbability = 0.25,
                         samplingStrategy = "uniform",
                         log = TRUE)

  return (lnl)

}

## ----echo = TRUE, eval=FALSE---------------------------------------------
#  set.seed(12345)
#  samplesUniform <- tess.mcmc(likelihoodFunction = likelihoodUniform,
#                              priors = priorsSampling,
#                              parameters = runif(2,0,1),
#                              logTransforms = c(TRUE,TRUE),
#                              delta = c(1,1),
#                              iterations = 10000,
#                              burnin = 1000,
#                              thinning = 10,
#                              adaptive = TRUE,
#                              verbose = TRUE)

## ----echo = FALSE, eval=TRUE---------------------------------------------
if ( file.exists("results/samplesUniform.rds") == FALSE ) {
set.seed(12345)
samplesUniform <- tess.mcmc(likelihoodFunction = likelihoodUniform,
                            priors = priorsSampling,
                            parameters = runif(2,0,1),
                            logTransforms = c(TRUE,TRUE),
                            delta = c(1,1),
                            iterations = 10000,
                            burnin = 1000,
                            thinning = 10,
                            adaptive = TRUE,
                            verbose = TRUE)
                            
saveRDS(samplesUniform, "results/samplesUniform.rds")
} else {
samplesUniform <- readRDS("results/samplesUniform.rds")
}

## ----label=mcmcPlotSamplingUniform, include=TRUE, fig.cap="Trace plots (left) and marginal posterior probability densities (right) for the diversification rate (speciation - extinction) and turnover rate (extinction) under \\emph{uniform} sampling from the MCMC simulation.", fig.pos='!ht'----
summary(samplesUniform)
plot(samplesUniform)

## ----echo=TRUE-----------------------------------------------------------
likelihoodDiversified <- function(params) {

  speciation <- params[1] + params[2]
  extinction <- params[2]

  lnl <- tess.likelihood(times.diversified,
                         lambda = speciation,
                         mu = extinction,
                         samplingProbability = 0.25,
                         samplingStrategy = "diversified",
                         log = TRUE)

  return (lnl)

}

## ----echo = TRUE, eval=FALSE---------------------------------------------
#  samplesDiversified <- tess.mcmc(likelihoodFunction = likelihoodDiversified,
#                              priors = priorsSampling,
#                              parameters = runif(2,0,1),
#                              logTransforms = c(TRUE,TRUE),
#                              delta = c(1,1),
#                              iterations = 10000,
#                              burnin = 1000,
#                              thinning = 10,
#                              adaptive = TRUE,
#                              verbose = TRUE)

## ----echo = FALSE, eval=TRUE---------------------------------------------
if ( file.exists("results/samplesDiversified.rds") == FALSE ) {

samplesDiversified <- tess.mcmc(likelihoodFunction = likelihoodDiversified,
                            priors = priorsSampling,
                            parameters = runif(2,0,1),
                            logTransforms = c(TRUE,TRUE),
                            delta = c(1,1),
                            iterations = 10000,
                            burnin = 1000,
                            thinning = 10,
                            adaptive = TRUE,
                            verbose = TRUE)
                            
saveRDS(samplesDiversified, "results/samplesDiversified.rds")
} else {
samplesDiversified <- readRDS("results/samplesDiversified.rds")
}

## ----label=mcmcPlotSamplingDiversified, include=TRUE, fig.cap="Trace plots (left) and marginal posterior probability densities (right) for the diversification rate (speciation - extinction) and turnover rate (extinction) under \\emph{diversified} sampling from the MCMC simulation.", fig.pos='!ht'----
summary(samplesDiversified)
plot(samplesDiversified)

## ----echo = TRUE, eval=FALSE---------------------------------------------
#  set.seed(12345)
#  marginalLikelihoodConstBD <- tess.steppingStoneSampling(
#                  likelihoodFunction = likelihoodConstBD,
#                  priors = priorsConstBD,
#                  parameters = runif(2,0,1),
#                  logTransforms = c(TRUE,TRUE),
#                  iterations = 1000,
#                  burnin = 100,
#                  K = 50)
#  
#  marginalLikelihoodDecrBD <- tess.steppingStoneSampling(
#                  likelihoodFunction = likelihoodDecrBD,
#                  priors = priorsDecrBD,
#                  parameters = runif(3,0,1),
#                  logTransforms = c(TRUE,TRUE,TRUE),
#                  iterations = 1000,
#                  burnin = 100,
#                  K = 50)
#  
#  marginalLikelihoodEpisodicBD <- tess.steppingStoneSampling(
#                  likelihoodFunction = likelihoodEpisodicBD,
#                  priors = priorsEpisodicBD,
#                  parameters = runif(4,0,1),
#                  logTransforms = c(TRUE,TRUE,TRUE,TRUE),
#                  iterations = 1000,
#                  burnin = 100,
#                  K = 50)
#  
#  marginalLikelihoodMassExtinctionBD <- tess.steppingStoneSampling(
#                  likelihoodFunction = likelihoodMassExtinctionBD,
#                  priors = priorsMassExtinctionBD,
#                  parameters = c(runif(2,0,1),max(times)*3/4),
#                  logTransforms = c(TRUE,TRUE,FALSE),
#                  iterations = 1000,
#                  burnin = 100,
#                  K = 50)
#  

## ----echo = FALSE, eval=TRUE---------------------------------------------
if ( file.exists("results/marginalLikelihoodConstBD.rds") == FALSE ) {

set.seed(12345)
marginalLikelihoodConstBD <- tess.steppingStoneSampling(
                likelihoodFunction = likelihoodConstBD,
                priors = priorsConstBD,
                parameters = runif(2,0,1),
                logTransforms = c(TRUE,TRUE),
                iterations = 1000,
                burnin = 100,
                K = 50)

marginalLikelihoodDecrBD <- tess.steppingStoneSampling(
                likelihoodFunction = likelihoodDecrBD,
                priors = priorsDecrBD,
                parameters = runif(3,0,1),
                logTransforms = c(TRUE,TRUE,TRUE),
                iterations = 1000,
                burnin = 100,
                K = 50)

marginalLikelihoodEpisodicBD <- tess.steppingStoneSampling(
                likelihoodFunction = likelihoodEpisodicBD,
                priors = priorsEpisodicBD,
                parameters = runif(4,0,1),
                logTransforms = c(TRUE,TRUE,TRUE,TRUE),
                iterations = 1000,
                burnin = 100,
                K = 50)

marginalLikelihoodMassExtinctionBD <- tess.steppingStoneSampling(
                likelihoodFunction = likelihoodMassExtinctionBD,
                priors = priorsMassExtinctionBD,
                parameters = c(runif(2,0,1),max(times)*3/4),
                logTransforms = c(TRUE,TRUE,FALSE),
                iterations = 1000,
                burnin = 100,
                K = 50)

saveRDS(marginalLikelihoodConstBD, "results/marginalLikelihoodConstBD.rds")
saveRDS(marginalLikelihoodDecrBD, "results/marginalLikelihoodDecrBD.rds")
saveRDS(marginalLikelihoodEpisodicBD, "results/marginalLikelihoodEpisodicBD.rds")
saveRDS(marginalLikelihoodMassExtinctionBD, "results/marginalLikelihoodMassExtinctionBD.rds")
} else {
marginalLikelihoodConstBD <- readRDS("results/marginalLikelihoodConstBD.rds")
marginalLikelihoodDecrBD <- readRDS("results/marginalLikelihoodDecrBD.rds")
marginalLikelihoodEpisodicBD <- readRDS("results/marginalLikelihoodEpisodicBD.rds")
marginalLikelihoodMassExtinctionBD <- readRDS("results/marginalLikelihoodMassExtinctionBD.rds")
}

## ----echo=TRUE-----------------------------------------------------------
# First, construct a vector of the marginal likelhoods named by the
# model to which they refer.
candidateModels <- c("ConstBD"=marginalLikelihoodConstBD,
                     "DecrBD"=marginalLikelihoodDecrBD,
                     "EpisodicBD"=marginalLikelihoodEpisodicBD,
                     "MassExtinctionBD"=marginalLikelihoodMassExtinctionBD)

# Make all possible combinations of the models.
marginalLikelihoodGrid <- expand.grid(M0=names(candidateModels),
                                      M1=names(candidateModels))

# Add a column that is the 2 ln BF for each pair of models.
marginalLikelihoodGrid$BF <- 2 * (candidateModels[marginalLikelihoodGrid$M0] -
                                  candidateModels[marginalLikelihoodGrid$M1])

# Sort the comparisons by their 2 ln BF in descending order.
marginalLikelihoodGrid <- marginalLikelihoodGrid[order(marginalLikelihoodGrid$BF,
                                                       decreasing=TRUE),]

marginalLikelihoodGrid

## ----echo = TRUE---------------------------------------------------------
tmrca <- max( times )

## ----echo = TRUE---------------------------------------------------------
# The simulation function
simConstBD <- function(params) {

  # Same model as above.
  speciation <- params[1] + params[2]
  extinction <- params[2]

  # We need trees with at least three tips for the
  # gamma-statistic.
  repeat {
    tree <- tess.sim.age(n = 1,
                         age = tmrca,
                         lambda = speciation,
                         mu = extinction,
                         samplingProbability = 1.0,
                         MRCA = TRUE)[[1]]
    if (tree$Nnode > 1) break
  }
  return (tree)
}

## ----echo = TRUE, eval=TRUE----------------------------------------------
# simulate trees from the posterior-predictive distribution
treesConstBD <- tess.PosteriorPrediction(simConstBD,samplesConstBD)

## ----echo = TRUE---------------------------------------------------------
# compute the number of species in each simulate tree
numTaxaConstBD <- c()
for (i in 1:length(treesConstBD)){
  numTaxaConstBD[i] <- treesConstBD[[i]]$Nnode + 1
}

## ----label = constPosteriorPredictiveNTaxa, echo = TRUE, include = TRUE, eval=FALSE----
#  # Compute the 95% posterior-predictive interval of the
#  # number of taxa.
#  numTaxaPPDI <- quantile(numTaxaConstBD,prob=c(0.025,0.975))
#  
#  # Plot the posterior-predictive distribution with the
#  # quantiles. Then, compare it to the observed number
#  # of species, x.
#  plot(density(numTaxaConstBD),main="Number of taxa",xlab="",
#      ylab="Posterior Predictive Density",lwd=2)
#  
#  abline(v=numTaxaPPDI,lty=2,col="gray",lwd=2)
#  
#  points(conifers$Nnode+1,0,pch="x")

## ----label = constPosteriorPredictiveLTT, echo = TRUE, include = TRUE, eval=FALSE----
#  ltt.plot(treesConstBD[[1]],backward=FALSE,col="gray",log="y",
#           ylim=c(1,max(numTaxaConstBD)),main="LTT-plot")
#  
#  for (i in 2:min(100,length(treesConstBD))) ltt.lines(treesConstBD[[i]],
#                             backward=FALSE, col="gray")
#  
#  ltt.lines(conifers,backward=FALSE,lwd=3)

## ----label = constPosteriorPredictiveGamma, echo = TRUE, include = TRUE, eval=FALSE----
#  
#  # Compute the observed gamma statistic.
#  observedGamma <- gammaStat(conifers)
#  
#  # Perform the posterior predictive test, and compute
#  # the 95% posterior predictive interval.
#  ppt <- tess.PosteriorPredictiveTest(treesConstBD,conifers,
#                                      gammaStat)
#  gammaPPDI <- quantile(ppt[[1]],prob=c(0.025,0.975))
#  
#  # Compare the observed statistic to the posterior
#  # predictive density.
#  plot(density(ppt[[1]]),main="Gamma Statistic",xlab="",
#                         ylab="Posterior Predictive Density",lwd=2)
#  abline(v=gammaPPDI,lty=2,col="gray",lwd=2)
#  points(observedGamma,0,pch="x")

## ----echo=FALSE, label=posteriorPredictiveTests, include=TRUE, fig.height=4, out.width="\\linewidth", fig.align="center", fig.cap="Assessing the absolute fit of the conifer tree to the constant-rate birth-death model using posterior-predictive simulation. (A) The posterior-predictive distribution for the number of species; the dashed gray lines indicate the 95\\% credible interval, and the `x' indicates the location of the observed species number. (B) LTT plots for the simulated trees (gray) and for the conifer study tree (black). (C) The posterior-predictive distribution for the gamma statistic; the dashed  gray lines indicate the 95\\% credible interval, and the `x' indicates the location of the value of the gamma statistic calculated for the conifer tree."----

par(mfrow=c(1,3),mar=c(5,4,3,0.1),las=1)

# Compute the 95% posterior predictive interval of the number of taxa.
numTaxaPosteriorPredictiveInterval <- quantile(numTaxaConstBD,prob=c(0.025,0.975))

# Plot the posterior predictive distribution with the quantiles.
# Then, compare it to the observed number of species, x.
plot(density(numTaxaConstBD),main=NA,xlab="",ylab="Posterior Predictive Density",lwd=2)
abline(v=numTaxaPosteriorPredictiveInterval,lty=2,col="gray",lwd=2)
points(conifers$Nnode+1,0,pch="x")
mtext("A", line = 1)

ltt.plot(treesConstBD[[1]],backward=FALSE,col="gray",log="y",ylim=c(1,max(numTaxaConstBD)),main=NA)
for (i in 2:min(100,length(treesConstBD))) ltt.lines(treesConstBD[[i]],backward=FALSE,col="gray")
ltt.lines(conifers,backward=FALSE,lwd=3)
mtext("B", line = 1)

# Compute the observed gamma statistic.
observedGamma <- gammaStat(conifers)

# Perform the posterior predictive test, and compute the 95% posterior predictive interval.
ppt <- tess.PosteriorPredictiveTest(treesConstBD,conifers,gammaStat)
gammaPosteriorPredictiveInterval <- quantile(ppt[[1]],prob=c(0.025,0.975))

# Compare the observed statistic to the posterior predictive density.
plot(density(ppt[[1]]),main=NA,xlab="",ylab="Posterior Predictive Density",lwd=2)
abline(v=gammaPosteriorPredictiveInterval,lty=2,col="gray",lwd=2)
points(observedGamma,0,pch="x")
mtext("C", line = 1)

## ----echo=TRUE, eval=TRUE------------------------------------------------
mean(ppt[[1]] >= observedGamma)

## ----echo = TRUE---------------------------------------------------------
numExpectedMassExtinctions <- 2

## ----echo=TRUE-----------------------------------------------------------
numExpectedRateChanges <- 2

## ----echo = TRUE---------------------------------------------------------
# Specify the mean and standard deviation of the lognormal
# prior on the speciation rate in real space
speciationPriorMu <- 0.2
speciationPriorSigma <- 0.5

# Specify the mean and standard deviation of the lognormal
# prior on the extinction rate in real space
extinctionPriorMu <- 0.15
extinctionPriorSigma <- 0.5

## ----echo = TRUE---------------------------------------------------------
# Transform the priors on the speciation rate into log space.
speciationRatePriorMean <- log((speciationPriorMu^2)
                           /sqrt(speciationPriorSigma^2+
                           speciationPriorMu^2))

speciationRatePriorStDev <- sqrt( log(1+speciationPriorSigma^2
                            /(speciationPriorMu^2)))

# Transform the priors on the extinction rate into log space.
extinctionRatePriorMean <- log((extinctionPriorMu^2)
                           /sqrt(extinctionPriorSigma^2+
                           extinctionPriorMu^2))

extinctionRatePriorStDev <- sqrt( log(1+extinctionPriorSigma^2
                            /(extinctionPriorMu^2)))

## ----echo=TRUE-----------------------------------------------------------
expectedSurvivalProbability <- 0.05

## ----echo=TRUE-----------------------------------------------------------
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 *
                                expectedSurvivalProbability /
                                (expectedSurvivalProbability - 1)

## ----echo = TRUE, label = priorSurvivalProbability, fig.cap="Our prior density on the survival probability of a mass-extinction event.",fig.height=3.5----
# Plot the density function of our beta distribution.
curve(dbeta(x,shape1=pMassExtinctionPriorShape1,
            shape2=pMassExtinctionPriorShape2),n=1001,
            xlab='survival probability',ylab='density',las=1)

# Plot the 95% prior interval on the survival probability.
abline(v = qbeta(c(0.025,0.975),shape1=pMassExtinctionPriorShape1,
                 shape2=pMassExtinctionPriorShape2),lty=2)

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  set.seed(12345)
#  tess.analysis(conifers,
#                empiricalHyperPriors = FALSE,
#                initialSpeciationRate = speciationPriorMu,
#                speciationRatePriorMean = speciationRatePriorMean,
#                speciationRatePriorStDev = speciationRatePriorStDev,
#                initialExtinctionRate = extinctionPriorMu,
#                extinctionRatePriorMean = extinctionRatePriorMean,
#                extinctionRatePriorStDev = extinctionRatePriorStDev,
#                samplingProbability = samplingFraction,
#                numExpectedRateChanges = numExpectedRateChanges,
#                numExpectedMassExtinctions = numExpectedMassExtinctions,
#                pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
#                pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
#                MAX_ITERATIONS = 10000,
#                dir = "tess_analysis")

## ----echo = FALSE, eval = TRUE-------------------------------------------
if ( file.exists("results/tess_analysis.rds") == FALSE ) {

set.seed(12345)
tess.analysis(conifers,
              empiricalHyperPriors = FALSE,
              initialSpeciationRate = speciationPriorMu,
              speciationRatePriorMean = speciationRatePriorMean,
              speciationRatePriorStDev = speciationRatePriorStDev,
              initialExtinctionRate = extinctionPriorMu,
              extinctionRatePriorMean = extinctionRatePriorMean,
              extinctionRatePriorStDev = extinctionRatePriorStDev,
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges,
              numExpectedMassExtinctions = numExpectedMassExtinctions,
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
              MAX_ITERATIONS = 10000,
              dir = "tess_analysis")
              
}

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  output <- tess.process.output("tess_analysis",
#                  numExpectedRateChanges = numExpectedRateChanges,
#                  numExpectedMassExtinctions = numExpectedMassExtinctions)
#  

## ----echo = FALSE, eval = TRUE-------------------------------------------
if ( file.exists("results/tess_analysis.rds") == FALSE ) {

output <- tess.process.output("tess_analysis",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)

saveRDS(output, "results/tess_analysis.rds")
} else {
output <- readRDS("results/tess_analysis.rds")
}

## ----echo = TRUE, eval = TRUE, label = plotTessAnalysis, fig.cap="Visualizing the results of a \\CoMET analysis when diversification hyperpriors are specified \\emph{a priori}."----

layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.output(output,
                 fig.types = c("speciation rates",
                               "speciation shift times",
                               "extinction rates",
                               "extinction shift times",
                               "mass extinction Bayes factors",
                               "mass extinction times"),
                 las=2)

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  set.seed(12345)
#  tess.analysis(conifers,
#                empiricalHyperPriors = TRUE,
#                samplingProbability = samplingFraction,
#                numExpectedRateChanges = numExpectedRateChanges,
#                numExpectedMassExtinctions = numExpectedMassExtinctions,
#                pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
#                pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
#                MAX_ITERATIONS = 10000,
#                dir = "comet_hyperpriors")
#  

## ----echo = FALSE, eval = TRUE-------------------------------------------
if ( file.exists("results/comet_hyperpriors.rds") == FALSE ) {

set.seed(12345)
tess.analysis(conifers,
              empiricalHyperPriors = TRUE,
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges,
              numExpectedMassExtinctions = numExpectedMassExtinctions,
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
              MAX_ITERATIONS = 10000,
              dir = "comet_hyperpriors")

}

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  output <- tess.process.output("comet_hyperpriors",
#                  numExpectedRateChanges = numExpectedRateChanges,
#                  numExpectedMassExtinctions = numExpectedMassExtinctions)
#  

## ----echo = FALSE, eval = TRUE-------------------------------------------
if ( file.exists("results/comet_hyperpriors.rds") == FALSE ) {

output <- tess.process.output("comet_hyperpriors",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)

saveRDS(output, "results/comet_hyperpriors.rds")
} else {
output <- readRDS("results/comet_hyperpriors.rds")
}

## ----echo = TRUE, eval = TRUE, label = plotTessAnalysisEmpirical, fig.cap="Visualizing the results of a \\CoMET analysis with empirically estimated diversification hyperpriors."----

layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.output(output,
                 fig.types = c("speciation rates",
                               "speciation shift times",
                               "extinction rates",
                               "extinction shift times",
                               "mass extinction Bayes factors",
                               "mass extinction times"),
                 las=2)

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  set.seed(12345)
#  tess.analysis(conifers,
#                empiricalHyperPriors = TRUE,
#                samplingProbability = samplingFraction,
#                estimateNumberRateChanges = FALSE,
#                numExpectedMassExtinctions = numExpectedMassExtinctions,
#                pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
#                pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
#                MAX_ITERATIONS = 10000,
#                dir = "comet_no_rateshifts")
#  
#  output <- tess.process.output("comet_no_rateshifts",
#                  numExpectedRateChanges = numExpectedRateChanges,
#                  numExpectedMassExtinctions = numExpectedMassExtinctions)

## ----echo = FALSE, eval = TRUE-------------------------------------------
if ( file.exists("results/comet_no_rateshifts.rds") == FALSE ) {

set.seed(12345)
tess.analysis(conifers,
              empiricalHyperPriors = TRUE,
              samplingProbability = samplingFraction,
              estimateNumberRateChanges = FALSE,
              numExpectedMassExtinctions = numExpectedMassExtinctions,
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
              MAX_ITERATIONS = 10000,
              dir = "comet_no_rateshifts")
              
output <- tess.process.output("comet_no_rateshifts",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)
                
saveRDS(output, "results/comet_no_rateshifts.rds")
} else {
output <- readRDS("results/comet_no_rateshifts.rds")

}

## ----echo = TRUE, eval = TRUE, label=plotTessAnalysisWithoutRateshifts, fig.cap = "Visualizing the results of a \\CoMET analysis with empirically estimated diversification hyperpriors \\emph{and without} diversification rate-shifts."----

layout.mat <- matrix(1:2,nrow=2,ncol=1)
layout(layout.mat)
tess.plot.output(output,
                 fig.types = c("mass extinction Bayes factors",
                               "mass extinction times"),
                 las=2)

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  set.seed(12345)
#  tess.analysis(conifers,
#                empiricalHyperPriors = TRUE,
#                samplingProbability = samplingFraction,
#                estimateNumberMassExtinctions = FALSE,
#                MAX_ITERATIONS = 10000,
#                dir = "comet_no_mass_extinctions")
#  
#  output <- tess.process.output("comet_no_mass_extinctions",
#                  numExpectedRateChanges = numExpectedRateChanges,
#                  numExpectedMassExtinctions = numExpectedMassExtinctions)

## ----echo = FALSE, eval = TRUE-------------------------------------------
if ( file.exists("results/comet_no_mass_extinctions.rds") == FALSE ) {

set.seed(12345)
tess.analysis(conifers,
              empiricalHyperPriors = TRUE,
              samplingProbability = samplingFraction,
              estimateNumberMassExtinctions = FALSE,
              MAX_ITERATIONS = 10000,
              dir = "comet_no_mass_extinctions")
              
output <- tess.process.output("comet_no_mass_extinctions",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)


saveRDS(output, "results/comet_no_mass_extinctions.rds")
} else {
output <- readRDS("results/comet_no_mass_extinctions.rds")
}

## ----echo = TRUE, eval = TRUE, label = plotTessAnalysisWithoutMassExtinctions, fig.cap="Visualizing the results of a \\CoMET analysis with empirically estimated diversification hyperpriors \\emph{and without} mass-extinction events."----
layout.mat <- matrix(1:4,nrow=2,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.output(output,
                 fig.types = c("speciation rates",
                               "speciation shift times",
                               "extinction rates",
                               "extinction shift times"),
                 las=2)

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
#  prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
#  my_priors <- c("diversification"=prior_delta,
#                     "turnover"=prior_tau)
#  my_likelihood <- function(params) {
#  
#    speciation <- params[1] + params[2]
#    extinction <- params[2]
#  
#    lnl <- tess.likelihood(times,
#                           lambda = speciation,
#                           mu = extinction,
#                           samplingProbability = 1.0,
#                           log = TRUE)
#  
#    return (lnl)
#  
#  }
#  
#  samples_run_1 <- tess.mcmc(likelihoodFunction = my_likelihood,
#                              priors = my_priors,
#                              parameters = runif(2,0,10),
#                              logTransforms = c(TRUE,TRUE),
#                              delta = c(1,1),
#                              iterations = 200,
#                              burnin = 0,
#                              thinning = 1,
#                              adaptive = TRUE,
#                              verbose = TRUE)
#  
#  samples_run_2 <- tess.mcmc(likelihoodFunction = my_likelihood,
#                              priors = my_priors,
#                              parameters = runif(2,0,10),
#                              logTransforms = c(TRUE,TRUE),
#                              delta = c(1,1),
#                              iterations = 200,
#                              burnin = 0,
#                              thinning = 1,
#                              adaptive = TRUE,
#                              verbose = TRUE)

## ----echo = FALSE, eval = TRUE-------------------------------------------
if ( file.exists("results/samples_run_1.rds") == FALSE ) {

prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
my_priors <- c("diversification"=prior_delta,
                   "turnover"=prior_tau)
my_likelihood <- function(params) {

  speciation <- params[1] + params[2]
  extinction <- params[2]

  lnl <- tess.likelihood(times,
                         lambda = speciation,
                         mu = extinction,
                         samplingProbability = 1.0,
                         log = TRUE)

  return (lnl)

}

samples_run_1 <- tess.mcmc(likelihoodFunction = my_likelihood,
                            priors = my_priors,
                            parameters = runif(2,0,10),
                            logTransforms = c(TRUE,TRUE),
                            delta = c(1,1),
                            iterations = 200,
                            burnin = 0,
                            thinning = 1,
                            adaptive = TRUE,
                            verbose = TRUE)

samples_run_2 <- tess.mcmc(likelihoodFunction = my_likelihood,
                            priors = my_priors,
                            parameters = runif(2,0,10),
                            logTransforms = c(TRUE,TRUE),
                            delta = c(1,1),
                            iterations = 200,
                            burnin = 0,
                            thinning = 1,
                            adaptive = TRUE,
                            verbose = TRUE)

saveRDS(samples_run_1, "results/samples_run_1.rds")
saveRDS(samples_run_2, "results/samples_run_2.rds")
} else {
samples_run_1 <- readRDS("results/samples_run_1.rds")
samples_run_2 <- readRDS("results/samples_run_2.rds")
}

## ----echo = TRUE---------------------------------------------------------
effectiveSize(samples_run_1)
effectiveSize(samples_run_2)

## ----echo = TRUE---------------------------------------------------------
geweke.diag(samples_run_1)
geweke.diag(samples_run_2)

## ----echo = TRUE---------------------------------------------------------
qnorm(0.05/2)

## ----echo = TRUE---------------------------------------------------------
qnorm(1-0.05/2)

## ----echo = TRUE---------------------------------------------------------
gelman.diag(x=list(run1 = samples_run_1,run2 = samples_run_2),
            confidence = 0.95, transform = FALSE,
            autoburnin = FALSE, multivariate=TRUE)$mpsrf

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  output <- tess.process.output("comet_hyperpriors",
#                  numExpectedRateChanges = numExpectedRateChanges,
#                  numExpectedMassExtinctions = numExpectedMassExtinctions)

## ----echo = FALSE, eval = TRUE-------------------------------------------
if ( file.exists("results/comet_hyperpriors.rds") == FALSE ) {

output <- tess.process.output("comet_hyperpriors",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)
                
saveRDS(output, "results/comet_hyperpriors.rds")
} else {
output <- readRDS("results/comet_hyperpriors.rds")
}

## ----echo = TRUE, eval = TRUE--------------------------------------------
# Compute the effective sample size and Geweke diagnostic for
# the number of speciation-rate shifts.
effectiveSize(output$numSpeciationCategories)
geweke.diag(output$numSpeciationCategories)

# Compute the effective sample size and Geweke diagnostic for
# the number of extinction-rate shifts.
effectiveSize(output$numExtinctionCategories)
geweke.diag(output$numExtinctionCategories)

# Compute the effective sample size and Geweke diagnostic for
# the number of mass-extinctionevents.
effectiveSize(output$numMassExtinctions)
geweke.diag(output$numMassExtinctions)

## ----echo = TRUE, label = tessAnalysisMCMCDiagnosis, eval = TRUE, fig.cap="Visualizing the single-chain MCMC diagnostics for a \\CoMET analysis with empirically estimated diversification hyperpriors. Blue bars/dots represent passed tests and red bars/dots mean failed tests (failed convergence)."----
layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.singlechain.diagnostics(output,
                   parameters = c("speciation rates",
                                  "extinction rates",
                                  "mass extinction times"),
                   las=2)

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  set.seed(12345)
#  posterior_directories <- paste("comet_posterior_",
#                                 1:4,sep="")
#  
#  for(dir in posterior_directories) {
#    tess.analysis(conifers,
#                  empiricalHyperPriors = TRUE,
#                  samplingProbability = samplingFraction,
#                  numExpectedRateChanges = numExpectedRateChanges,
#                  numExpectedMassExtinctions = numExpectedMassExtinctions,
#                  pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
#                  pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
#                  MAX_ITERATIONS = 10000,
#                  dir = dir)
#  }

## ----echo = FALSE, results = "hide", eval = TRUE-------------------------
if ( file.exists("results/comet_posterior_1.rds") == FALSE ) {

set.seed(12345)
posterior_directories <- paste("comet_posterior_",
                               1:4,sep="")

for(dir in posterior_directories) {
  tess.analysis(conifers,
                empiricalHyperPriors = TRUE,
                samplingProbability = samplingFraction,
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions,
                pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
                pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
                MAX_ITERATIONS = 10000,
                dir = dir)
                
}

output_1 <- tess.process.output("comet_posterior_1",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)

output_2 <- tess.process.output("comet_posterior_2",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)

output_3 <- tess.process.output("comet_posterior_3",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)

output_4 <- tess.process.output("comet_posterior_4",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)



saveRDS(output_1, "results/comet_posterior_1.rds")
saveRDS(output_2, "results/comet_posterior_2.rds")
saveRDS(output_3, "results/comet_posterior_3.rds")
saveRDS(output_4, "results/comet_posterior_4.rds")

} else {

output_1 <- readRDS("results/comet_posterior_1.rds")
output_2 <- readRDS("results/comet_posterior_2.rds")
output_3 <- readRDS("results/comet_posterior_3.rds")
output_4 <- readRDS("results/comet_posterior_4.rds")

}

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  output_1 <- tess.process.output("comet_posterior_1",
#                  numExpectedRateChanges = numExpectedRateChanges,
#                  numExpectedMassExtinctions = numExpectedMassExtinctions)
#  
#  output_2 <- tess.process.output("comet_posterior_2",
#                  numExpectedRateChanges = numExpectedRateChanges,
#                  numExpectedMassExtinctions = numExpectedMassExtinctions)
#  
#  output_3 <- tess.process.output("comet_posterior_3",
#                  numExpectedRateChanges = numExpectedRateChanges,
#                  numExpectedMassExtinctions = numExpectedMassExtinctions)
#  
#  output_4 <- tess.process.output("comet_posterior_4",
#                  numExpectedRateChanges = numExpectedRateChanges,
#                  numExpectedMassExtinctions = numExpectedMassExtinctions)

## ----echo = TRUE, eval = TRUE, fig.cap="Visualizing the multiple-chain MCMC diagnostics for a \\CoMET analysis with empirically estimated diversification hyperpriors. Blue dots represent passed tests and red dots mean failed tests (failed convergence)."----
output_list <- list(output_1,output_2,output_3,output_4)

layout.mat <- matrix(1:3,nrow=3,ncol=1,byrow=TRUE)
layout(layout.mat)
tess.plot.multichain.diagnostics(output_list,
                                 parameters = c("speciation rates",
                                                "extinction rates",
                                                "mass extinction times"),
                                 las=2)

## ----echo = TRUE---------------------------------------------------------
prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
my_priors <- c("diversification"=prior_delta,
                   "turnover"=prior_tau)
my_likelihood <- function(params) {

  speciation <- params[1] + params[2]
  extinction <- params[2]

  lnl <- tess.likelihood(times,
                         lambda = speciation,
                         mu = extinction,
                         samplingProbability = 1.0,
                         log = TRUE)

  return (lnl)

}

## ----echo = TRUE, eval= FALSE--------------------------------------------
#  samples_run_small <- tess.mcmc(likelihoodFunction = my_likelihood,
#                              priors = my_priors,
#                              parameters = runif(2,0,10),
#                              logTransforms = c(TRUE,TRUE),
#                              delta = c(0.1,0.02),
#                              iterations = 1000,
#                              burnin = 200,
#                              thinning = 1,
#                              adaptive = FALSE,
#                              verbose = TRUE)
#  

## ----echo = FALSE, eval = TRUE-------------------------------------------
if ( file.exists("results/samples_run_small.rds") == FALSE ) {

samples_run_small <- tess.mcmc(likelihoodFunction = my_likelihood,
                            priors = my_priors,
                            parameters = runif(2,0,10),
                            logTransforms = c(TRUE,TRUE),
                            delta = c(0.1,0.02),
                            iterations = 1000,
                            burnin = 200,
                            thinning = 1,
                            adaptive = FALSE,
                            verbose = TRUE)

saveRDS(samples_run_small, "results/samples_run_small.rds")
} else {
samples_run_small <- readRDS("results/samples_run_small.rds")
}

## ----label=mcmcNoTuningPlotSmall, include=TRUE, fig.cap="Trace plots (left) and marginal posterior probability densities (right) for the diversification rate (top) and turnover rate (bottom) from the MCMC simulation under the constant-rate birth-death process.", fig.pos='!ht'----
plot(samples_run_small)

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  samples_run_good <- tess.mcmc(likelihoodFunction = my_likelihood,
#                              priors = my_priors,
#                              parameters = runif(2,0,10),
#                              logTransforms = c(TRUE,TRUE),
#                              delta = c(1,0.2),
#                              iterations = 1000,
#                              burnin = 200,
#                              thinning = 1,
#                              adaptive = FALSE,
#                              verbose = TRUE)
#  

## ----echo = FALSE, eval = TRUE-------------------------------------------
if ( file.exists("results/samples_run_good.rds") == FALSE ) {

samples_run_good <- tess.mcmc(likelihoodFunction = my_likelihood,
                            priors = my_priors,
                            parameters = runif(2,0,10),
                            logTransforms = c(TRUE,TRUE),
                            delta = c(1,0.2),
                            iterations = 1000,
                            burnin = 200,
                            thinning = 1,
                            adaptive = FALSE,
                            verbose = TRUE)

saveRDS(samples_run_good, "results/samples_run_good.rds")
} else {
samples_run_good <- readRDS("results/samples_run_good.rds")
}

## ----label=mcmcNoTuningPlotGood, include=TRUE, fig.cap="Trace plots (left) and marginal posterior probability densities (right) for the diversification rate (top) and turnover rate (bottom) from the MCMC simulation under the constant-rate birth-death process.", fig.pos='!ht'----
plot(samples_run_good)

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  samples_run_large <- tess.mcmc(likelihoodFunction = my_likelihood,
#                              priors = my_priors,
#                              parameters = runif(2,0,10),
#                              logTransforms = c(TRUE,TRUE),
#                              delta = c(10,2),
#                              iterations = 1000,
#                              burnin = 200,
#                              thinning = 1,
#                              adaptive = FALSE,
#                              verbose = TRUE)

## ----echo = FALSE, eval = TRUE-------------------------------------------
if ( file.exists("results/samples_run_large.rds") == FALSE ) {

samples_run_large <- tess.mcmc(likelihoodFunction = my_likelihood,
                            priors = my_priors,
                            parameters = runif(2,0,10),
                            logTransforms = c(TRUE,TRUE),
                            delta = c(10,2),
                            iterations = 1000,
                            burnin = 200,
                            thinning = 1,
                            adaptive = FALSE,
                            verbose = TRUE)

saveRDS(samples_run_large, "results/samples_run_large.rds")
} else {
samples_run_large <- readRDS("results/samples_run_large.rds")
}

## ----label=mcmcNoTuningPlotLarge, include=TRUE, fig.cap="Trace plots (left) and marginal posterior probability densities (right) for the diversification rate (top) and turnover rate (bottom) from the MCMC simulation under the constant-rate birth-death process.", fig.pos='!ht'----
plot(samples_run_large)

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  samples_run_auto_tuned <- tess.mcmc(likelihoodFunction = my_likelihood,
#                              priors = my_priors,
#                              parameters = runif(2,0,10),
#                              logTransforms = c(TRUE,TRUE),
#                              delta = c(1,1),
#                              iterations = 1000,
#                              burnin = 200,
#                              thinning = 1,
#                              adaptive = TRUE,
#                              verbose = TRUE)

## ----echo = FALSE, eval = TRUE-------------------------------------------
if ( file.exists("results/samples_run_auto_tuned.rds") == FALSE ) {

samples_run_auto_tuned <- tess.mcmc(likelihoodFunction = my_likelihood,
                            priors = my_priors,
                            parameters = runif(2,0,10),
                            logTransforms = c(TRUE,TRUE),
                            delta = c(1,1),
                            iterations = 1000,
                            burnin = 200,
                            thinning = 1,
                            adaptive = TRUE,
                            verbose = TRUE)
                            
saveRDS(samples_run_auto_tuned, "results/samples_run_auto_tuned.rds")
} else {
samples_run_auto_tuned <- readRDS("results/samples_run_auto_tuned.rds")
}

## ----label=mcmcTuningPlot, include=TRUE, fig.cap="Trace plots (left) and marginal posterior probability densities (right) for the diversification rate (top) and turnover rate (bottom) from the MCMC simulation under the constant-rate birth-death process.", fig.pos='!ht'----
plot(samples_run_auto_tuned)

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  set.seed(12345)
#  
#  tess.analysis(conifers,
#                empiricalHyperPriors = TRUE,
#                samplingProbability = samplingFraction,
#                numExpectedRateChanges = numExpectedRateChanges,
#                numExpectedMassExtinctions = numExpectedMassExtinctions,
#                pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
#                pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
#                MAX_ITERATIONS = 100000000,
#                MAX_TIME = 24*60*60,
#                MIN_ESS = 500,
#                dir = "comet_auto_stop")
#  
#  
#  output_auto_stop <- tess.process.output("comet_auto_stop",
#            numExpectedRateChanges = numExpectedRateChanges,
#            numExpectedMassExtinctions = numExpectedMassExtinctions)

## ----echo = FALSE, results = "hide", eval = TRUE-------------------------
if ( file.exists("results/output_auto_stop.rds") == FALSE ) {

set.seed(12345)

tess.analysis(conifers,
              empiricalHyperPriors = TRUE,
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges,
              numExpectedMassExtinctions = numExpectedMassExtinctions,
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
              MAX_ITERATIONS = 100000000,
              MAX_TIME = 24*60*60,
              MIN_ESS = 500,
              dir = "comet_auto_stop")
              
              
output_auto_stop <- tess.process.output("comet_auto_stop",
          numExpectedRateChanges = numExpectedRateChanges,
          numExpectedMassExtinctions = numExpectedMassExtinctions)
              
saveRDS(output_auto_stop, "results/output_auto_stop.rds")
} else {
output_auto_stop <- readRDS("results/output_auto_stop.rds")
}

## ----echo = TRUE, eval = TRUE, label = cometAutoStopPlot-----------------

# Compute the effective sample size and Geweke diagnostic for
# the number of speciation-rate shifts.
effectiveSize(output_auto_stop$numSpeciationCategories)
geweke.diag(output_auto_stop$numSpeciationCategories)

# Compute the effective sample size and Geweke diagnostic for
# the number of extinction-rate shifts.
effectiveSize(output_auto_stop$numExtinctionCategories)
geweke.diag(output_auto_stop$numExtinctionCategories)

# Compute the effective sample size and Geweke diagnostic for
# the number of mass-extinction events.
effectiveSize(output_auto_stop$numMassExtinctions)
geweke.diag(output_auto_stop$numMassExtinctions)

layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.singlechain.diagnostics(output_auto_stop,
                   parameters = c("speciation rates",
                                  "extinction rates",
                                  "mass extinction times"),
                   las=2)

