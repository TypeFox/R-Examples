rm(list = ls())
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
library(PAWL)

# fix seed
set.seed(17)
# define a mixture model
ncomponents <- 4
MP <- list(
  ncomponents = ncomponents,
  componentweights = rep(0.25, ncomponents), 
  componentmeans = c(-3, 0, 3, 6),
  componentvariances = rep(0.55^2, ncomponents))
mixture <- createMixtureTarget(mixturesize = 100, 
                               ncomponents = ncomponents, 
                               mixtureparameters = MP)


schematictargets <- data.frame(cbind(c(rep(-3, 3), rep(0, 3), rep(3, 3), rep(6, 3)), 
                                     c(0, 3, 6, -3, 3, 6, -3, 0, 6, -3, 0, 3)))
g <- ggplot(schematictargets, aes(x = X1, y = X2))
g <- g + geom_point(size = 20, colour = "red")
g <- g + geom_point(size = 15, colour = "white")
g <- g + geom_point(size = 10, colour = "red")
g <- g + geom_point(size = 5, colour = "white")
g <- g + xlim(-5, 8) + ylim(-5, 8)
g <- g + xlab(expression(mu[1])) + ylab(expression(mu[2]))
# where are the modes ?
print(g)

# save this mixture model (because we are going to alter it slightly during the script)
originalmixture <- mixture

###### launch PAWL

# set the algorithmic parameters
# number of chains
N <- 10
# number of iterations
T <- 10^4
# burnin (set it to 1 if you don't want any burnin)
burnin <- 10^3
# we are going to store all the chains in this example
mcmcparameters <- tuningparameters(nchains = N, niterations = T, storeall = TRUE,
                                   computemean = TRUE, computemeanburnin = burnin)

cat("PAWL starting...\n")
ptm <- proc.time()
preexpresults <- preexplorationAMH(target = mixture, nchains = mcmcparameters@nchains, 
                                   niterations = 1000, verbose = TRUE)

binrange <- preexpresults$SuggestedRange
cat("bin range:", binrange, "\n")
rinit <- function(size) preexpresults$finalchains
mixture@rinit <- rinit
ncuts <- 20
getPos <- function(points, logdensity) - logdensity 
energybinning <- binning(position = getPos, name = "energy",
                         binrange = binrange, ncuts = ncuts,
                         useLearningRate = TRUE, autobinning = TRUE,
                         alongenergy = TRUE)
pawlresults <- pawl(target = mixture, binning = energybinning, 
                    AP = mcmcparameters, verbose = TRUE)
runtime <- proc.time() - ptm

cat("# FH met:", length(pawlresults$FHtimes), "\n")
#print(apply(pawlresults$meanchains, 2, mean))
pawlchains <- ConvertResults(pawlresults, verbose = FALSE)
locations <- energybinning@getLocations(pawlresults$finalbins, - pawlchains$logdens)
finallogtheta <- pawlresults$logtheta[dim(pawlresults$logtheta)[1],]
finaltheta <- exp(finallogtheta - max(finallogtheta))
finaltheta <- finaltheta / sum(finaltheta)
pawlchains$importanceweights <- finaltheta[locations]
pawlchains <- subset(pawlchains, iterations > mcmcparameters@computemeanburnin)
pawlchains <- subset(pawlchains, select = c("X5", "X6", "indexchain", 
                                            "iterations", "logdens", "importanceweights"))
names(pawlchains) <- c("Mu1", "Mu2", "indexchain", "iterations", "logdens", "importanceweights")
totalnpoints <- dim(pawlchains)[1]
pawlchains$index <- 1:totalnpoints
maxnumberpoints <- min(totalnpoints, 50000)
pawlchains <- subset(pawlchains, index %in% sample(1:totalnpoints, maxnumberpoints, replace = FALSE))
g <- ggplot(pawlchains, aes(x = Mu1, y = Mu2))
g <- g + geom_point(aes(alpha = importanceweights))
g <- g + xlab(expression(mu[1])) + ylab(expression(mu[2]))
g <- g + theme(legend.position = "none")
g <- g + scale_alpha(range=c(0.05, 0.1))
g <- g + xlim(-5, 8) + ylim(-5, 8)
# plot (mu_1, mu_2)
print(g)

###### launch parallel adaptive MH for comparison
# since we modified mixture@rinit, we go back to the original target
# before doing anything else
mixture <- originalmixture
#  let's use the same algorithmic parameters

ptm <- proc.time()
amhresults <- adaptiveMH(mixture, mcmcparameters, verbose = TRUE)
runtime <- proc.time() - ptm
amhchains <- ConvertResults(amhresults, verbose = FALSE)
amhchains <- subset(amhchains, select = c("X5", "X6", "indexchain", "iterations", "logdens"))
names(amhchains) <- c("Mu1", "Mu2", "indexchain", "iterations", "logdens")
amhchains <- subset(amhchains, iterations > burnin)
totalnpoints <- dim(amhchains)[1]
amhchains$index <- 1:totalnpoints
maxnumberpoints <- min(totalnpoints, 500000)
amhchains <- subset(amhchains, index %in% sample(1:totalnpoints, maxnumberpoints, replace = FALSE))
g <- ggplot(amhchains, aes(x = Mu1, y = Mu2))
g <- g + geom_point(alpha = 1/20)
g <- g + xlab(expression(mu[1])) + ylab(expression(mu[2]))
g <- g + theme(legend.position = "none")
g <- g + xlim(-5, 8) + ylim(-5, 8)
print(g)

###### launch SMC
mixture <- originalmixture
smcparameters <- smcparameters(nparticles=5000, 
                               temperatures = seq(from = 0.0001, to = 1, length.out= 100),
                               nmoves = 5, ESSthreshold = 0.5, movetype = "randomwalk",
                               movescale = 0.1)
ptm <- proc.time()
smcresults <- smc(mixture, smcparameters, verbose = TRUE)
runtime <- proc.time() - ptm
x <- smcresults$particles
w <- normalizeweight(smcresults$weights)
particles <- data.frame(cbind(x[,5:8], w))
names(particles) <- c("Mu1", "Mu2", "Mu3", "Mu4", "weights")
g <- ggplot(particles, aes(x = Mu1, y = Mu2))
g <- g + geom_point(aes(alpha = weights))
g <- g + xlab(expression(mu[1])) + ylab(expression(mu[2]))
g <- g + theme(legend.position = "none")
g <- g + scale_alpha(range=c(0.05, 0.1))
g <- g + xlim(-5, 8) + ylim(-5, 8)
print(g)

