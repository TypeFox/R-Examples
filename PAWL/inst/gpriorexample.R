# remove all objects
rm(list = ls())
# setting a seed for the RNG
set.seed(17)
# try to detach the package if it was already loaded
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
# load the package
library(PAWL)
library(fields)

data(pollution)
Yvector = Pollution[,1]
Xmatrix = Pollution[,2:16]

### the g factor is set to exp(20)
g <- exp(20)
targetdimension <- dim(Xmatrix)[2];
temperature = 1;

targetparameters <- list(g = g, y = Yvector, 
                         X = Xmatrix, n = length(Yvector), 
                         ysquared = t(Yvector) %*% Yvector,
                         targetdimension = targetdimension,
                         temperature = temperature)


##### initialization 
rInitDistribution <- function(size){
  matrix(sample(c(1,0), targetdimension*size, replace=T), nr=size)
}

lpostw <- function(gammaVS, tp){  
    # Xt1 is the matrix of the selected covariates (selected gammaVS == 1)
    Xt1 <- as.matrix(tp$X[, gammaVS==1])
    # creating n x n matrix denoted by P1
    if (sum(gammaVS) != 0){
        #P1 <- Xt1 %*% solve(t(Xt1) %*% Xt1) %*% t(Xt1)
        # using crossprod, it's a bit faster
        P1 <- Xt1 %*% tcrossprod(x = solve(crossprod(x = Xt1, y = Xt1)), y = Xt1)
    } else {
        P1 <- matrix(0, tp$n, tp$n)
    }
    # computing the result, in two parts since the formula is quite long
    posteriorvalue <- -(sum(gammaVS) + 1) / 2 * log(tp$g + 1) 
    posteriorvalue <- posteriorvalue - tp$n / 2 * log(tp$ysquared - tp$g / (tp$g+1) * t(tp$y) %*% P1 %*% tp$y) 
    return(posteriorvalue)
}

### target density function of the algorithm
logdensity <- function(gammaVS, targetparameters){
    if (is.vector(gammaVS)){
        logpost <- lpostw(gammaVS, targetparameters)
    } else {
        logpost <- rep(0, dim(gammaVS)[1])
        for(j in 1:dim(gammaVS)[1]){
            logpost[j] <- lpostw(gammaVS[j,], targetparameters)
        }
    }
    return(targetparameters$temperature*logpost)
}
    

gpriortarget <- target(name = "gprior", type = "discrete", dimension = targetdimension,
                         rinit = rInitDistribution, logdensity = logdensity,
                         parameters = targetparameters)
print(gpriortarget)

proposalparam <- list(targetdimension = targetdimension)

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

#preexpresults <- preexplorationAMH(gpriortarget, 5, 100000, proposalinstance)
#binrange <- preexpresults$SuggestedRange
binrange <- c(377,450)
#rinit <- function(size)
#  preexpresults$finalchains
#gpriortarget@rinit <- rinit

getLogEnergy <- function(points, logdensity) -logdensity
densitybinning <- binning(position = getLogEnergy,
                          name = "minus log target density",
                          binrange = binrange,
                          ncuts = 20,
                          autobinning = FALSE,
                          alongenergy = FALSE)
                          
nchains = 1; T = 80000;
mhparameters1 <- tuningparameters(nchains = nchains, niterations = T, 
                                 saveeverynth = T/500, computemean = TRUE)
print(mhparameters1)
now = proc.time()
pawlresults1 <- pawl(gpriortarget, densitybinning, mhparameters1, proposalinstance)
proc.time() - now

nchains = 10; T = 25000;
mhparameters10 <- tuningparameters(nchains = nchains, niterations = T, 
                                 saveeverynth = T/500, computemean = TRUE)
print(mhparameters10)
now = proc.time()
pawlresults10 <- pawl(gpriortarget, densitybinning, mhparameters10, proposalinstance)
proc.time() - now

nchains = 100; T = 3500;
mhparameters100 <- tuningparameters(nchains = nchains, niterations = T, 
                                 saveeverynth = T/500, computemean = TRUE)
print(mhparameters100)
now = proc.time()
pawlresults100 <- pawl(gpriortarget, densitybinning, mhparameters100, proposalinstance)
proc.time() - now
                                   
nchains = 100; T = 3500
amhparameters <- tuningparameters(nchains = nchains, niterations = T, 
                                 saveeverynth = T/500, computemean = TRUE)
amhresults1 <- adaptiveMH(gpriortarget, amhparameters, proposalinstance)

gpriortarget@parameters$temperature = 0.1;
amhresults0_1 <- adaptiveMH(gpriortarget, amhparameters, proposalinstance)

gpriortarget@parameters$temperature = 0.01;
amhresults0_01 <- adaptiveMH(gpriortarget, amhparameters, proposalinstance)


                         
# And now, a ton of code to produce pretty plots...
                         
                         
meltIt = function(particles, nbparticles, name){
  p_mat = as.matrix(apply(particles, c(1,2),mean))
	p_mean = apply(p_mat, 1, mean)
	p_upper = apply(p_mat, 1, quantile, probs=.95)
	p_lower = apply(p_mat, 1, quantile, probs=.05)
	return(cbind(melt(rep(name,length(p_mean))), melt(1:length(p_mean)), melt(p_mean), melt(p_lower)$value, melt(p_upper)$value))
}

meltedpawl1 = meltIt (pawlresults1$allchains, 1, "Wang-Landau")
meltedpawl10 = meltIt (pawlresults10$allchains, 10, "Wang-Landau")
meltedpawl100 = meltIt (pawlresults100$allchains, 100, "Wang-Landau")

meltedamh1 = meltIt (amhresults1$allchains, 100, "Metropolis-Hastings, Temp = 1")
meltedamh0_1 = meltIt (amhresults0_1$allchains, 100, "Metropolis-Hastings, Temp = 10")
meltedamh0_01 = meltIt (amhresults0_01$allchains, 100, "Metropolis-Hastings, Temp = 100")

melted <- rbind(meltedpawl100, meltedamh1, meltedamh0_1, meltedamh0_01)
names(melted) <- c("N","Iteration","Mean", "Lower", "Upper")
melted$N = as.factor(melted$N)

pdf(file = "gPriorExploration.pdf", width = 10, height = 10)
p <- ggplot(melted, aes(x = as.numeric(Iteration)*7, y = Mean, ymin = Lower, ymax = Upper, colour = N, fill = N)) 
p <- p + geom_ribbon(alpha = .4, aes(colour = NULL)) + geom_line(size = 2) + ylab("Model Saturation") + xlab("Iteration")
p <- p + facet_wrap( ~ N, ncol=2) + theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12, angle = 90) ) + ylim(0,.75)
print(p)
dev.off()

                         
                         
## Convergence plot

#Grab true values (obtained from enumeration)
trueProportions = c(9.344054e-01, 3.510045e-02, 2.797744e-02, 2.376642e-03, 1.313087e-04, 8.407919e-06, 3.532125e-07, 1.292003e-08, 3.676946e-10, 1.469434e-11, 3.640064e-13, 1.251748e-14, 2.586274e-16, 6.659019e-18, 1.269542e-19, 2.706594e-21, 6.558651e-23, 9.459380e-25, 2.310697e-26, 2.387800e-28, 6.134736e-30)

normalizeTheta <- function(logtheta){
  logtheta - log(sum(exp(logtheta)))
}

meltThetaValues = function(logtheta, name, thinningFactor){
  for(i in 1:dim(logtheta)[1]){
    logtheta[i,] = normalizeTheta(logtheta[i,])
  }
  RMSE = sqrt(rowMeans((logtheta- matrix(log(trueProportions), byrow=T, nr=dim(logtheta)[1], nc=dim(logtheta)[2]))^2))[seq(1,dim(logtheta)[1],by=thinningFactor)]
  meltedConvergence = melt(logtheta[seq(1,dim(logtheta)[1],by=thinningFactor),])
  meltedConvergence <- cbind(name, meltedConvergence)
  names(meltedConvergence) = c("N", "Iteration","Bin","LogTheta")
  meltedConvergence$Iteration = meltedConvergence$Iteration*thinningFactor
  return(list(meltedConvergence = meltedConvergence, RMSE=RMSE))
}

cnvg1 = meltThetaValues(pawlresults1$logtheta, "N = 1", 80000/500)
cnvg10 = meltThetaValues(pawlresults10$logtheta, "N = 10", 25000/500)
cnvg100 = meltThetaValues(pawlresults100$logtheta, "N = 100", 3500/500)

meltedConvergence = rbind(cnvg1$meltedConvergence, cnvg10$meltedConvergence, cnvg100$meltedConvergence)

RMSE = melt(cbind(cnvg1$RMSE, cnvg10$RMSE, cnvg100$RMSE))
names(RMSE) = c("Iteration","N","RMSE")
RMSE$N[RMSE$N == 1] <- "N=1"
RMSE$N[RMSE$N == 2] <- "N=10"
RMSE$N[RMSE$N == 3] <- "N=100"


pdf(file = "gPriorConvergence.pdf", width = 10, height = 4)
p <- ggplot(meltedConvergence, aes(Iteration, LogTheta, group = Bin))
p <- p + geom_line(aes(color = Bin)) + geom_hline(aes(yintercept=log(trueProportions), color = 1:length(trueProportions)), linetype="dotted")
p <- p + xlab("Iteration") + ylab(expression(paste("Log(",theta,")"))) + theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12, angle = 90) )
p <- p + facet_grid(. ~ N, scales="free") + ylim(-75,5)
print(p)
dev.off()