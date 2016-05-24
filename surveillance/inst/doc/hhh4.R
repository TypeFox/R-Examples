### R code from vignette source 'hhh4.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: setup
###################################################
library("surveillance")  
options(width=75)

## create directory for plots
dir.create("plots", showWarnings=FALSE)

######################################################
## Do we need to compute or can we just fetch results?
######################################################
compute <- !file.exists("hhh4-cache.RData")
message("Doing computations: ", compute)
if(!compute) load("hhh4-cache.RData")


###################################################
### code chunk number 2: loadInfluMen
###################################################
# load data
data("influMen")
# convert to sts class and print basic information about the time series
print(fluMen <- disProg2sts(influMen))


###################################################
### code chunk number 3: getMen
###################################################
meningo <- fluMen[, "meningococcus"]
dim(meningo)


###################################################
### code chunk number 4: plotfluMen
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(fluMen, type = observed ~ time | unit, # type of plot (default)
             same.scale = FALSE,            # unit-specific ylim?
             col = "grey")                  # color of bars


###################################################
### code chunk number 5: readInFlu
###################################################
# read in observed number of cases
flu.counts <- as.matrix(read.table(system.file("extdata/counts_flu_BYBW.txt", 
                                               package = "surveillance"),
                                   check.names = FALSE))


###################################################
### code chunk number 6: nhoodByBw
###################################################
getOption("SweaveHooks")[["fig"]]()
# read in 0/1 adjacency matrix (1 if regions share a common border)
nhood <- as.matrix(read.table(system.file("extdata/neighbourhood_BYBW.txt",
                                          package = "surveillance"),
                              check.names = FALSE))
library("Matrix")
print(image(Matrix(nhood)))


###################################################
### code chunk number 7: fluAsSTS
###################################################
# read in population fractions
popfracs <- read.table(system.file("extdata/population_2001-12-31_BYBW.txt",
                                   package = "surveillance"),
                       header = TRUE)$popFrac
# create sts object
flu <- sts(flu.counts, start = c(2001, 1), frequency = 52,
           population = popfracs, neighbourhood = nhood)


###################################################
### code chunk number 8: plot-flu-ByBw
###################################################
getOption("SweaveHooks")[["fig"]]()
data("fluBYBW")
plot(fluBYBW[year(fluBYBW) == 2001, ], # select year 2001
     type = observed ~ unit,           # total counts by region
     population = fluBYBW@map$X31_12_01 / 100000) # per 100000 inhabitants
grid::grid.text("Incidence [per 100'000 inhabitants]", x = 0.5, y = 0.02)


###################################################
### code chunk number 9: hhh4.Rnw:271-276
###################################################
# consistency check
local({
    fluBYBW@map <- flu@map
    stopifnot(all.equal(fluBYBW, flu))
})


###################################################
### code chunk number 10: measles2w
###################################################
data("measlesDE")
measles2w <- aggregate(measlesDE, nfreq = 26)


###################################################
### code chunk number 11: plot-measles
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(measles2w, type = observed ~ time,  # aggregate counts over all units
     main = "Bi-weekly number of measles cases in Germany")


###################################################
### code chunk number 12: hhh4 (eval = FALSE)
###################################################
## hhh4(sts, control)


###################################################
### code chunk number 13: controlObj (eval = FALSE)
###################################################
## control = list(
##     ar = list(f = ~ -1,              # formula for log(lambda_it)
##               offset = 1),           # optional multiplicative offset
##     ne = list(f = ~ -1,              # formula for log(phi_it)
##               offset = 1,            # optional multiplicative offset
##               weights = neighbourhood(stsObj) == 1),  # (w_ji) matrix
##     end = list(f = ~ 1,              # formula for log(nu_it)
##                offset = 1),          # optional multiplicative offset e_it
##     family = "Poisson",              # Poisson or NegBin model
##     subset = 2:nrow(stsObj),         # subset of observations to be used 
##     optimizer = list(stop = list(tol = 1e-5, niter = 100), # stop rules
##                      regression = list(method = "nlminb"), # for penLogLik
##                      variance = list(method = "nlminb")),  # for marLogLik
##     verbose = FALSE,                 # level of progress reporting
##     start = list(fixed = NULL,       # list with initial values for fixed,
##                  random = NULL,      # random, and
##                  sd.corr = NULL),    # variance parameters
##     data = list(t = epoch(stsObj)-1),# named list of covariates 
##     keep.terms = FALSE               # whether to keep the model terms
## )


###################################################
### code chunk number 14: fitMeningo0
###################################################
# specify a formula object for the endemic component
( f_S1 <- addSeason2formula(f = ~ 1, S = 1, period = 52) )
# fit the Poisson model
result0 <- hhh4(meningo, control = list(end = list(f = f_S1),
                                        family = "Poisson"))
summary(result0)


###################################################
### code chunk number 15: fitMeningo1
###################################################
result1 <- update(result0, family = "NegBin1")


###################################################
### code chunk number 16: hhh4.Rnw:501-502
###################################################
AIC(result0, result1)


###################################################
### code chunk number 17: fitMeningo2
###################################################
# fit an autoregressive model
result2 <- update(result1, ar = list(f = ~ 1))


###################################################
### code chunk number 18: hhh4.Rnw:515-519
###################################################
coef(result2, se = TRUE,    # also return standard errors
     amplitudeShift = TRUE, # transform sine/cosine coefficients
                            # to amplitude/shift parameters
     idx2Exp = TRUE)        # exponentiate remaining parameters


###################################################
### code chunk number 19: plot_result2
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(result2)


###################################################
### code chunk number 20: neighbourhood_fluMen
###################################################
# no "transmission" from meningococcus to influenza
neighbourhood(fluMen)["meningococcus","influenza"] <- 0
neighbourhood(fluMen)


###################################################
### code chunk number 21: fitFluMen
###################################################
# create formula for endemic component
f.end <- addSeason2formula(f = ~ -1 + fe(1, unitSpecific = TRUE),
                                           # disease-specific intercepts
                           S = c(3, 1),    # S = 3 for flu, S = 1 for men
                           period = 52)
# specify model
m <- list(ar = list(f = ~ -1 + fe(1, unitSpecific = TRUE)),
          ne = list(f = ~ 1,  # phi, only relevant for meningococcus due to
                    weights = neighbourhood(fluMen)),   # the weight matrix
          end = list(f = f.end),
          family = "NegBinM") # disease-specific overdispersion
# fit model
result <- hhh4(fluMen, control = m)
summary(result, idx2Exp=1:3)


###################################################
### code chunk number 22: plot-fit_men
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(result, units = 1:2, legend = 2, legend.args = list(
     legend = c("influenza-driven", "autoregressive", "endemic")))


###################################################
### code chunk number 23: ri (eval = FALSE)
###################################################
## f.end <- ~ -1 + ri(type = "iid", corr = "all")


###################################################
### code chunk number 24: modelFluBYBW
###################################################
# endemic component: iid random effects, linear trend, S=3 seasonal terms
f.end <- addSeason2formula(f = ~ -1 + ri(type="iid", corr="all") +
                               I((t-208)/100),
                           S = 3, period = 52)
# model specification
model.B2 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ -1 + ri(type="iid", corr="all"),
                           weights = neighbourhood(fluBYBW),
                           normalize = TRUE),  # all(rowSums(weights) == 1)
                 end = list(f = f.end, offset = population(fluBYBW)),
                 family = "NegBin1", verbose = TRUE,
                 optimizer = list(variance = list(method = "Nelder-Mead")))
# default start values for random effects are sampled from a normal
set.seed(42)


###################################################
### code chunk number 25: computeFluBYBW
###################################################
if(compute){
  result.B2 <- hhh4(fluBYBW, model.B2)
  s.B2 <- summary(result.B2, maxEV = TRUE, idx2Exp = 1:3)
  
  #pred.B2 <- oneStepAhead(result.B2, tp = nrow(fluBYBW) - 2*52)
  predfinal.B2 <- oneStepAhead(result.B2, tp = nrow(fluBYBW) - 2*52,
                               type = "final")
  meanSc.B2 <- colMeans(scores(predfinal.B2))
  
  save(s.B2, meanSc.B2, file="hhh4-cache.RData")
}


###################################################
### code chunk number 26: fitFluBYBW (eval = FALSE)
###################################################
## # fit the model (takes about 35 seconds)
## result.B2 <- hhh4(fluBYBW, model.B2)
## summary(result.B2, maxEV = TRUE, idx2Exp = 1:3)


###################################################
### code chunk number 27: hhh4.Rnw:666-667
###################################################
s.B2


###################################################
### code chunk number 28: oneStepAhead_rolling (eval = FALSE)
###################################################
## pred.B2 <- oneStepAhead(result.B2, tp = nrow(fluBYBW) - 2*52)


###################################################
### code chunk number 29: oneStepAhead_fake (eval = FALSE)
###################################################
## predfinal.B2 <- oneStepAhead(result.B2, tp = nrow(fluBYBW) - 2*52,
##                              type = "final")


###################################################
### code chunk number 30: scores (eval = FALSE)
###################################################
## colMeans(scores(predfinal.B2, which = c("logs", "rps")))


###################################################
### code chunk number 31: hhh4.Rnw:699-700
###################################################
meanSc.B2[c("logs", "rps")]


###################################################
### code chunk number 32: createVacc
###################################################
data(MMRcoverageDE)
cardVac1 <- MMRcoverageDE[1:16,3:4]

adjustVac <- function(cardVac, p=0.5,nrow=1){
  card <- cardVac[,1]
  vac <- cardVac[,2]
  vacAdj <- vac*card + p*vac*(1-card)
  return(matrix(vacAdj,nrow=nrow, ncol=length(vacAdj), byrow=TRUE))
}
vac0 <- 1-adjustVac(cardVac1,p=0.5,nrow=measles2w@freq*3)
colnames(vac0) <- colnames(measles2w)


###################################################
### code chunk number 33: hhh4.Rnw:746-747
###################################################
vac0[1:2, 1:6]


###################################################
### code chunk number 34: fitMeasles
###################################################
# endemic component: Intercept + sine/cosine terms
f.end <- addSeason2formula(f = ~ 1, S = 1, period = 26)
# autoregressive component: Intercept + vaccination coverage information
model.A0 <- list(ar = list(f = ~ 1 + logVac0),
                 end = list(f = f.end, offset = population(measles2w)),
                 data = list(t = epoch(measles2w), logVac0 = log(vac0)))
# fit the model
result.A0 <- hhh4(measles2w, model.A0)
summary(result.A0, amplitudeShift = TRUE)


