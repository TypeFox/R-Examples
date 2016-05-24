# remove all objects
#graphics.off()
#rm(list = ls())
# try to detach the package if it was already loaded
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
# load the package
library(PAWL)
# starting points for MCMC algorithms
rinit <- function(size) rnorm(size)
# target log density function: a gaussian distribution N(mean = 2, sd = 3)
parameters <- list(mean = 2, sd = 3)
logdensity <- function(x, parameters) dnorm(x, parameters$mean, parameters$sd, log = TRUE)
# creating the target object
gaussiantarget <- target(name = "gaussian", dimension = 1,
                         rinit = rinit, logdensity = logdensity,
                         parameters = parameters)
# setting a seed for the RNG
set.seed(17)

#######
## Adaptive Metropolis-Hastings
#######
#mhparameters <- tuningparameters(nchains = 10, niterations = 1000, storeall = TRUE)
#amhresults <- adaptiveMH(gaussiantarget, mhparameters)
#range(amhresults$alllogtarget)
## check that it's working
#PlotHist(results = amhresults, component = 1)
#curve(dnorm(x, mean = gaussiantarget@parameters$mean,
#            sd = gaussiantarget@parameters$sd), add = TRUE, lwd = 2, col = "red")


######
# Parallel Adaptive Wang-Landau
######
N <- 10
T <- 2000

# here we disable the adaptive proposal to highlight the automatic binning
# mechanism
proposal <- createAdaptiveRandomWalkProposal(N, gaussiantarget@dimension,
                                             adaptiveproposal = FALSE)

pawlparameters <- tuningparameters(nchains = N, niterations = T, storeall = TRUE)
print(pawlparameters)

#########
# PAWL where we bin along the state space
getPos <- function(points, logdensity) points
## we further specify some parameters, like the bins,
## the desired frequency in each bin
ncuts <- 2
positionbinning <- binning(position = getPos,
                            name = "position",
                            binrange = c(-20, -3),
                            ncuts = ncuts,
                            useLearningRate = TRUE,
                            autobinning = TRUE)

print(positionbinning)
pawlresults <- pawl(gaussiantarget, binning = positionbinning, AP = pawlparameters,
                    proposal = proposal)
# histogram of the binned coordinate
PlotHistBin(pawlresults, positionbinning)
#########
## PAWL where we bin along the energy (= - log density)
getPos <- function(points, logdensity) - logdensity 

# we can get the range using pre exploratory mcmc
preexpresults <- preexplorationAMH(gaussiantarget, N, 1000)
binrange <- preexpresults$SuggestedRange
# or specify it ourselves
#binrange <- c(2.1, 20)
# we further specify some parameters, like the bins,
# the desired frequency in each bin...
ncuts <- 2
energybinning <- binning(position = getPos,
                            name = "energy",
                            binrange = binrange,
                            ncuts = ncuts,
                            useLearningRate = TRUE,
                            autobinning = TRUE,
                            alongenergy = TRUE)


print(energybinning)
# launching the algorithm...
pawlresults <- pawl(gaussiantarget, binning = energybinning, AP = pawlparameters,
                    proposal = proposal)
# histogram of the binned coordinate
#PlotHistBin(pawlresults, energybinning)
##
#########

#########
## results
#
getFrequencies(pawlresults, energybinning)
# Histogram of the chains
#PlotHist(pawlresults, 1)
# Trace plot of the log theta penalties
#print(PlotLogTheta(pawlresults))
# trace plot of all the variables
# (here there is only one variable)
#print(PlotAllVar(pawlresults))







