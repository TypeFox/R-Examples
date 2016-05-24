library(ggplot2)
theme_update(
axis.title.x = element_text(size=25),
axis.title.y = element_text(size=25, angle = 90),
axis.text.x = element_text(size=25),
axis.text.y = element_text(size=25),
strip.text.x = element_text(size=25),
strip.text.y = element_text(size=25),
plot.title = element_text(size=25),
legend.text = element_text(size=25),
legend.title = element_text(size=25),
strip.background = element_rect(fill = "whitesmoke"))


rm(list = ls())
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
library(PAWL)

set.seed(17)
trimodal <- createTrimodalTarget()
N <- 2
Tprelim <- 500
preexp = preexplorationAMH(target = trimodal, nchains = N, niterations = Tprelim)
print("Suggesting this energy range:")
print(preexp$SuggestedRange)

T <- 2500
getLogEnergy <- function(points, logdensity) -logdensity
densitybinning <- binning(position = getLogEnergy,
                            name = "minus log target density",
                            binrange = preexp$SuggestedRange,
                            ncuts = 2,
                            useFH = TRUE,
                            autobinning = TRUE,
                            smoothbinning = TRUE,
                            splitThreshold = 0.25)

print(densitybinning)
pawlparameters <- tuningparameters(nchains = N, niterations = T, storeall = TRUE)
print(pawlparameters)

### Launching the algorithm...
#rinit <- function(size)
#  preexp$finalchains
#trimodal@rinit <- rinit
pawlresults <- pawl(trimodal, binning = densitybinning, AP = pawlparameters)
getFrequencies(pawlresults, densitybinning)

chains <- ConvertResults(pawlresults)

# 2D density plot of the components
T <- max(chains$iterations)
burnin <- min(1000, T / 10)
subchains <- subset(chains, iterations > burnin)
totalnpoints <- dim(subchains)[1]
subchains$index <- 1:totalnpoints
maxnumberpoints <- 50000
subchains <- subset(subchains, index > totalnpoints - maxnumberpoints)
g <- ggplot(subchains, aes(x = X1, y = X2))
g <- g + stat_bin2d() + geom_density2d()
g <- g + theme(legend.position = "none")
g <- g + xlab(expression(X[1])) + ylab(expression(X[2]))
#pdf(file = "Trimodal2Ddensity.pdf")
#print(g)
#dev.off()
# cloud of points with colour representing the density values
g <- ggplot(data = subchains, aes(x = X1, y = X2))
g <- g + geom_point(aes(alpha = logdens))  
g <- g + xlab(expression(X[1])) + ylab(expression(X[2]))
g <- g + theme(legend.position = "none")
g <- g + xlim(-15, 15) + ylim(-15, 15)
#ggsave(g, file = "TrimodalCloud.png")

# trace plot of log theta
st <- pawlresults$splitTimes
T <- length(pawlresults$acceptrates)
st <- c(0, st, T)
library(foreach)
library(reshape)
library(ggplot2)
df <- foreach (i= 1:(length(st) - 1), .combine = rbind) %do% {
    substart <- st[i] + 1
    substop  <- st[i+1] + 1
    sublogtheta <- data.frame(pawlresults$logthetahistory[[i]])
    #thetaDF <- exp(sublogtheta) / apply(exp(sublogtheta), 1, sum)
    thetaDF <- sublogtheta
    names(thetaDF) <- paste("theta", seq(1, pawlresults$nbins[i]))
    thetaDF$iterations <- substart:substop
    mdata <- melt(thetaDF, id = c("iterations"))
    names(mdata) <- c("iterations", "estimator", "value")
    mdata
}
# trace plot of log theta around the first split
g <- ggplot(subset(df, iterations < st[length(st)-1] + 200), aes(x = iterations, y = value, colour = estimator))
g <- ggplot(df, aes(x = iterations, y = value, colour = estimator))
g <- g + geom_line() 
g <- g + geom_vline(xintercept = pawlresults$splitTimes, linetype = 1)
g <- g + theme(legend.position = "none")
#g <- g + xlim(0, st[length(st)-1] + 200)
#pdf(file = "TrimodalLogThetasSplit.pdf", )
#pdf(file = "TrimodalLogThetasSplit.pdf", width = 21, height = 7)
#print(g)
#dev.off()

### We can get precise estimates
## of the true thetas, which are equal (in bin i) to:
##  psi_i / phi_i
## (renormalized)
## where psi_i is the integral of the target over bin i
## and phi_i is the desired frequency of the bin
proposedvalues <- trimodal@generate(10^6, trimodal@parameters)
logtde <- trimodal@logdensity(proposedvalues, trimodal@parameters)
proposedvalues <- cbind(proposedvalues, logtde)
BINS <- pawlresults$finalbins
locations <- densitybinning@getLocations(BINS, -proposedvalues[,"logtde"])
truethetas <- tabulate(locations, nbins = length(BINS))
truethetas <- truethetas / pawlresults$finaldesiredfreq
truethetas <- truethetas / sum(truethetas)

tabulate(densitybinning@getLocations(BINS, pawlresults$allreaction))

df <- foreach (i= 1:(length(st) - 1), .combine = rbind) %do% {
    substart <- st[i] + 1
    substop  <- st[i+1] + 1
    sublogtheta <- data.frame(pawlresults$logthetahistory[[i]])
    thetaDF <- exp(sublogtheta) / apply(exp(sublogtheta), 1, sum)
    names(thetaDF) <- paste("theta", seq(1, pawlresults$nbins[i]))
    thetaDF$iterations <- substart:substop
    mdata <- melt(thetaDF, id = c("iterations"))
    names(mdata) <- c("iterations", "estimator", "value")
    mdata
}
# trace plot of log theta around the first split
df$i <- 1:(dim(df)[1])
maxnumberpoints <- 10000
iterstep <- floor(dim(df)[1] / maxnumberpoints) + 1

# trace plot of log theta between the last
# bin split and the final iteration
g <- ggplot(subset(df, iterations %% iterstep == 0), aes(x = iterations, y = value, colour = estimator))
g <- g + geom_line() + scale_y_log10()
g <- g + geom_hline(yintercept = truethetas, linetype = 3)
g <- g + theme(legend.position = "none")
g <- g + xlim(st[length(st) - 1], T)
#pdf(file = "TrimodalLogThetasStable.pdf")
#print(g)
#dev.off()

# histogram of the energy values
Xnames <- grep("X", names(chains), value = TRUE)
positions <- data.frame(densitybinning@position(chains[,Xnames], chains$logdens))
names(positions) <- c("energy")
npoints <- dim(positions)[1]
positions$index <- 1:npoints
maxnumberpoints <- 500000
g <- ggplot(data = subset(positions, index > npoints - maxnumberpoints), aes(x = energy))
g <- g + geom_histogram(binwidth = 0.025, aes(y = ..density..))
g <- g + geom_vline(xintercept = densitybinning@bins[-1], size = 2)
g <- g + geom_vline(xintercept = pawlresults$finalbins[-1], linetype = 2, size = 2)
#pdf(file = "TrimodalHistogramBins.pdf", width = 21, height = 7)
#print(g)
#dev.off()


#### And now the adaptive MCMC with the same number of target density
#### evaluations. Since we don't use a preliminary exploration here,
#### the chains are run for N + Nprelim iterations.
##
mhparameters <- tuningparameters(nchains = N, niterations = T + Tprelim,
                                 storeall = TRUE) 
##
##### launching the algorithm...
amhresults <- adaptiveMH(trimodal, mhparameters)
#
#PlotAllVar(amhresults)
#pdf(file = "Marginal.pdf")
#par(mfrow = c(2, 1))
#PlotHist(amhresults, 1)
#PlotHist(pawlresults, 1)
#dev.off()
#pdf(file = "Trimodal2DdensityAMH.pdf")
#print(PlotDensComp1vsComp2(amhresults, "X1", "X2"))
#dev.off()

# 2D density plot of the components
amhchains <- ConvertResults(amhresults)
T <- max(amhchains$iterations)
burnin <- min(1000, T / 10)
subchains <- subset(amhchains, iterations > burnin)
totalnpoints <- dim(subchains)[1]
subchains$index <- 1:totalnpoints
maxnumberpoints <- 50000
subchains <- subset(subchains, index > totalnpoints - maxnumberpoints)
# cloud of points with colour representing the density values
g <- ggplot(data = subchains, aes(x = X1, y = X2))
g <- g + geom_point(aes(alpha = logdens))  
g <- g + xlab(expression(X[1])) + ylab(expression(X[2]))
g <- g + theme(legend.position = "none")
g <- g + xlim(-15, 15) + ylim(-15, 15)
#ggsave(g, file = "TrimodalCloudAMH.png")

