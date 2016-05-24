# remove all objects
rm(list = ls())
# setting a seed for the RNG
set.seed(17)
# try to detach the package if it was already loaded
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
# load the package
library(PAWL)
library(fields)
data(icefloe)

imgsize <- dim(IceFloe)[1]
targetdimension <- imgsize * imgsize
targetparameters <- list(a = 1, b = 0.8, imgmatrix = IceFloe, 
                         targetdimension = targetdimension,
                         imagesize = imgsize)

##### initialization 
rInitDistribution <- function(size){
  initpoints <- matrix(nrow = size, ncol = targetdimension)
  for(i in 1:size){
    initpoints[i,] <- c(IceFloe)
  }
  initpoints
}

logdensity <- function(X, parameters){
    out <- .Call("IsingCounter", chains = X, dataimg = parameters$imgmatrix,
                 imagesize = parameters$imagesize)
    countData <- out$countData
    countSimilar <- out$countSimilarities
    return(parameters$a * countData + parameters$b * (countSimilar / 2))
}
logdensityupdate <- function(chains, parameters, updateparam){
    res <- .Call("IsingUpdate", chains = chains, dataimg = parameters$imgmatrix,
                imagesize = parameters$imagesize, flippedindex = updateparam)
    updateLikelihood <- 2*(1/2 - res$equalToData) * parameters$a
    updatePrior <- res$SimilarityChange * parameters$b / 2
    return(updateLikelihood + updatePrior)
}
isingtarget <- target(name = "ising", type = "discrete", dimension = targetdimension,
                         rinit = rInitDistribution, logdensity = logdensity,
                         parameters = targetparameters, 
                         logdensityupdate = logdensityupdate)
print(isingtarget)

proposalparam <- list(imagesize = imgsize, targetdimension = targetdimension)

rproposal <- function(states, proposalparam){
    nchains <- dim(states)[1]
    index_to_flip <- sample.int(n = proposalparam$targetdimension, 
                                size = nchains, replace = TRUE)
    states[nchains * (index_to_flip - 1) + 1:nchains] <- !(states[nchains *
                                                           (index_to_flip - 1) + 1:nchains])
  return(list(states = states, others = index_to_flip))
}
# function to compute the density of the proposal kernel
# (necessary to compute the acceptance rate)
dproposal <- function(states, ys, proposalparam){
  return(rep(0, dim(states)[1]))
}

proposalinstance <- proposal(rproposal = rproposal, 
                             dproposal = dproposal,
                             proposalparam = proposalparam)

######
# Adaptive Metropolis-Hastings
######
#mhparameters <- tuningparameters(nchains = 10, niterations = 2500, 
#                                 storeall = FALSE)
N <- 10
T <- 10^6
mhparameters <- tuningparameters(nchains = N, niterations = T, 
                                 saveeverynth = T / 100, computemean = TRUE) 
print(mhparameters)

amh <- adaptiveMH(isingtarget, mhparameters, proposalinstance)

preexpresults <- preexplorationAMH(isingtarget, N, 20000, proposalinstance)
binrange <- preexpresults$SuggestedRange
rinit <- function(size)
  preexpresults$finalchains
isingtarget@rinit <- rinit

getLogEnergy <- function(points, logdensity) -logdensity
densitybinning <- binning(position = getLogEnergy,
                          name = "minus log target density",
                          binrange = binrange,
                          ncuts = 20,
                          autobinning = TRUE,
                          alongenergy = TRUE)
#
#
#Rprof(tmp <- tempfile())
pawlresults <- pawl(isingtarget, densitybinning, mhparameters, proposalinstance)
#Rprof()
# display profiling results
#print(summaryRprof(tmp))
#unlink(tmp)
# Plots showing exploration over time


library(ggplot2)
grid = expand.grid(1:40,1:40)


melted1 = rbind(cbind(rep(1,1600), rep(3,1600), grid[,1], grid[,2], 1-amh$allchains[60,2,]),
cbind(rep(1,1600), rep(4,1600), grid[,1], grid[,2], 1-amh$allchains[70,2,]),
cbind(rep(1,1600), rep(5,1600), grid[,1], grid[,2], 1-amh$allchains[80,2,]),
cbind(rep(1,1600), rep(6,1600), grid[,1], grid[,2], 1-amh$allchains[90,2,]),
cbind(rep(1,1600), rep(7,1600), grid[,1], grid[,2], 1-amh$allchains[100,2,]) )

melted2 = rbind(cbind(rep(2,1600), rep(3,1600), grid[,1], grid[,2], 1-pawlresults$allchains[60,2,]),
cbind(rep(2,1600), rep(4,1600), grid[,1], grid[,2], 1-pawlresults$allchains[70,2,]),
cbind(rep(2,1600), rep(5,1600), grid[,1], grid[,2], 1-pawlresults$allchains[80,2,]),
cbind(rep(2,1600), rep(6,1600), grid[,1], grid[,2], 1-pawlresults$allchains[90,2,]),
cbind(rep(2,1600), rep(7,1600), grid[,1], grid[,2], 1-pawlresults$allchains[100,2,]) )

melted = data.frame(rbind(melted1, melted2))
names(melted) = c("Group", "Sample", "X1", "X2", "Pixel")
melted$Sample[melted$Sample==3] = "Iteration 600,000"
melted$Sample[melted$Sample==4] = "Iteration 700,000"
melted$Sample[melted$Sample==5] = "Iteration 800,000"
melted$Sample[melted$Sample==6] = "Iteration 900,000"
melted$Sample[melted$Sample==7] = "Iteration 1,000,000"


melted$Group[melted$Group==1] = "Metropolis-Hastings"
melted$Group[melted$Group==2] = "Wang-Landau"
melted$Group = factor(melted$Group)

melted$Pixel[melted$Pixel==0] = "Off"
melted$Pixel[melted$Pixel==1] = "On"
melted$Pixel = factor(melted$Pixel)
melted$Pixel = factor(melted$Pixel, levels = rev(levels(melted$Pixel)))


pdf(file="AllExplore.pdf", width=12, height=5)
p <- ggplot(melted, aes(X1, X2, fill = Pixel)) 
p <- p + geom_tile() + facet_grid(Group ~ Sample)
p <- p + xlab(expression(X[1])) + ylab(expression(X[2])) + scale_fill_grey()
print(p)
dev.off()
           
           
meltedMean = rbind(cbind(rep(1,1600), grid[,1], grid[,2], apply(amh$allchains[60:100,,], 3, mean)), cbind(rep(2,1600), grid[,1], grid[,2], apply(pawlresults$allchains[60:100,,], 3, mean)))

meltedMean = data.frame(meltedMean)
names(meltedMean) = c("Group", "X1", "X2", "Pixel")

meltedMean$Group[meltedMean$Group==1] = "Metropolis-Hastings"
meltedMean$Group[meltedMean$Group==2] = "Wang-Landau"
meltedMean$Group = factor(meltedMean$Group)
      
#meltedMean$Pixel[meltedMean$Pixel==0] = "Off"
#meltedMean$Pixel[meltedMean$Pixel==1] = "On"
meltedMean$Pixel = factor(round(2*meltedMean$Pixel,digits=1)/2)
meltedMean$Pixel = factor(meltedMean$Pixel, levels = rev(levels(meltedMean$Pixel)))
meltedMean$Pixel = (as.numeric(meltedMean$Pixel))^(1/2) / (20)^(1/2)

pdf(file="MeanExplore.pdf", width=12, height=6)
p <- ggplot(meltedMean, aes(X1, X2, fill = Pixel)) 
p <- p + geom_tile() + facet_grid(. ~ Group) + xlab(expression(X[1])) + ylab(expression(X[2])) + scale_fill_gradient(low = "white", high = "black")
print(p)
dev.off()



