#rm(list = ls())
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
library(PAWL)

### toy example: the state space is made of three states
# target distribution:
targetpdf <- c(0.5, 0.4, 0.1)
# target log probability:
parameters <- list(logpi = log(targetpdf))
logdensity <- function(x, parameters){
    parameters$logpi[x]
}
# we have to specify a proposal mechanism
# for the Metropolis-Hastings kernel
# since the (default) continuous
# gaussian random walk cannot be used
transitionmatrix <- t(matrix(c(0.5, 0.5, 0.0,
                               0.3, 0.5, 0.2,
                               0.0, 0.6, 0.4), ncol = 3))

proposalparam <- list(transitionmatrix = transitionmatrix, card = 3)
# function that generates proposal:
rproposal <- function(states, proposalparam){
    for (index in 1:length(states)){
        states[index] <- sample(x = 1:proposalparam$card, 
                                size = 1, prob = proposalparam$transitionmatrix[states[index],])
    }
  return(list(states = states))
}
# function to compute the density of the proposal kernel
# (necessary to compute the acceptance rate)
dproposal <- function(states, ys, proposalparam){
    for (index in 1:(length(states))){
        states[index] <- log(transitionmatrix[states[index], ys[index]])
    }
  return(states)
}

proposalinstance <- proposal(rproposal = rproposal, 
                             dproposal = dproposal, 
                             proposalparam = proposalparam)

# function to draw starting points for the MCMC algorithms:
rinit <- function(size) return(rep(1, size))
# define the target
discretetarget <- target(name = "discrete toy example", dimension = 1, type = "discrete",
                         rinit = rinit, logdensity = logdensity, parameters = parameters)
# specify Metropolis-Hastings tuning parameters:
mhparameters <- tuningparameters(nchains = 1, niterations = 10000, storeall = TRUE)
# Rprof(tmp <- tempfile())
amhresults <- adaptiveMH(discretetarget, mhparameters, proposalinstance)
# Rprof()
# print(summaryRprof(tmp))
# unlink(tmp)
chains <- ConvertResults(amhresults)

cat("AMH: target probabilities:", targetpdf, "\n")
amhcount <- tabulate(chains$X1, nbins = 3)
cat("AMH: obtained frequencies:", amhcount / sum(amhcount), "\n")

# we bin such that states 1 and 2 are in bin 1, and state 3 is in bin 2
getPos <- function(points, logdensity) 2 - (points <= 2)
# we further specify some parameters, like the bins,
# the desired frequency in each bin...
positionbinning <- binning(position = getPos,
                            name = "position",
                            bins = c(1, 2),
                            desiredfreq = c(0.8, 0.2),
                            useLearningRate = FALSE)
pawlresults <- pawl(discretetarget, binning = positionbinning, AP = mhparameters, proposalinstance)
pawlchains <- ConvertResults(pawlresults)
cat("PAWL: desired frequencies:", positionbinning@desiredfreq, "\n")
pawlcount <- tabulate(getPos(pawlchains$X1, pawlchains$logdens), nbins = 2)
cat("PAWL: obtained frequencies:", pawlcount / sum(pawlcount), "\n")
# show the trace plot of log theta:
#PlotLogTheta(pawlresults)

counts <- tabulate(pawlchains$X1, nbins = 3)
counts <- counts[1:2]
counts <- counts / sum(counts)
cat("PAWL: obtained proportions inside bin 1:", counts, "\n")
cat("PAWL: compared to:", targetpdf[1:2] / sum(targetpdf[1:2]), "\n")
