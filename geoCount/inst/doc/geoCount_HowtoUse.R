

####################################################
### The code for the tutorial "How to Use geoCount"
### Author: Liang JIng (ljing918@gmail.com)
### edited on Dec 09 2013 to adapt geoCount 1.131209
####################################################

##################
### load packages
##################
require(geoCount)
require(coda)
require(reldist)
require(distrEx)
require(snowfall)
require(rlecuyer)
# require(multicore)
# {multicore} currently is not available in Windows!

################################################
################################################
### CHAPTER 2: Data Simulation and Visualization
################################################
################################################

#######################
### Simulate Locations
#######################
## grid locations
loc <- locGrid(1, 2, 10, 5)
plot(loc, xlab="x", ylab="y")

## circular locations
loc2 <- rbind( locCircle(1, 60), locCircle(0.667, 40), 
        locCircle(0.333, 20),  locCircle(0, 1) )
plot(loc2, xlab="x", ylab="y")

## squared locations
loc3 <- rbind( locSquad(1, 15), locSquad(0.667, 10), 
        locSquad(0.333, 5),  c(0,0) )
plot(loc3, xlab="x", ylab="y")

## scale locations to fit in a unit square
loc4 <- unifLoc(loc2, length=1)
plot(loc4, xlab="x", ylab="y")

#######################
### Simulate Data
#######################
loc <- locGrid(1, 1, 10, 10)
dat <- simData(loc = loc, L = 0, 
        X = NULL, beta = 4, cov.par = c(1, 0.2, 1), 
        rho.family = "rhoPowerExp", Y.family = "Poisson")
plotData(dat$data, loc, xlab="x", ylab="y")

#######################
### Data Visualization
#######################
loc <- rbind(locCircle(1, 60), 
               locCircle(0.667, 40), 
               locCircle(0.333, 20)  
               )
dat <- simData(loc, cov.par = c(1, 0.2, 1))
Y <- dat$data
  plotData(Y[1:60], loc[1:60, ], Y[61:100], loc[61:100, ], 
           Y[101:120], loc[101:120, ], pch = 16,
           xlab="x", ylab="y"
          )
## plot with boundaries
data(TexasCounty_boundary)
plotData(bdry = TexasCounty.boundary, xlab = "Longitude", ylab = "Latitude")

########################
### Integrated Data Sets
########################
## Rongelap data
data(Rongelap)
str(Rongelap)
plotData(bdry = Rongelap$borders, Y = Rongelap$data, loc = Rongelap$coords)

## Weed data
data(Weed)
str(Weed)
plotData(Weed[,3], Weed[,1:2], Weed[,4], Weed[,1:2], 
         xlab="East (meter)", ylab="North (meter)")

## Earthquakes data
data(Earthquake)
str(Earthquakes)
range(Earthquakes$Mag)

## TexasCounty data
data(TexasCounty_population)
str(TexasCounty.population)
data(TexasCounty_center)
str(TexasCounty.center)
data(TexasCounty_boundary)
length(TexasCounty.boundary)
plotTexas(TexasCounty.boundary) 

################################################
################################################
### CHAPTER 3: Estimation and Prediction
################################################
################################################

########################
### Environment Setting
########################
input <- MCMCinput( run = 10000, run.S = 10,
          rho.family = "rhoPowerExp",
          Y.family = "Poisson", ifkappa=0,
          scales=c(0.5, 1.5, 0.9, 0.6, 0.5),
          phi.bound=c(0.005, 1),
          initials=list(c(-1, 2, 1), 1, 0.1, 1) )

########################
### Run MCMC Algorithms
########################
## For example of analyzing Weed data,
  data(Weed)
  input.Weed <- MCMCinput( run=2000, run.S=1, rho.family="rhoPowerExp", 
          Y.family = "Poisson", ifkappa=0,
          scales=c(0.2, 3.5, 0.9, 0.6, 0.5), 
          phi.bound=c(0.5, 300), 
          initials=list(c(1), 1, 0.1, 1) )
  res <- runMCMC(Y=Weed[,3], L=0, loc=Weed[,1:2], X=NULL, MCMCinput=input.Weed )

## Include the covariates in the model
  input2.Weed <- MCMCinput( run=1000, run.S=1, rho.family="rhoPowerExp", 
          Y.family = "Poisson", ifkappa=0,
          scales=c(0.5, 0.00005, 0.9, 0.9, 0.5), 
          phi.bound=c(0.5, 300), 
          initials=list(c(4, 0, 0), 1, 0.1, 1) )
  res2 <- runMCMC(Y=Weed[,3], L=0, loc=Weed[,1:2], X=Weed[,1:2], 
                  MCMCinput=input2.Weed )

############################
### Generate Parallel Chains
############################
## For example, the parallel version of analyzing Weed data

# CAUSION: runMCMC.multiChain FUNCTION IS REMOVED since geoCount 1.131208 
# BECAUSE: any package depending on {multicore} is not accepted by CRAN any more.
# parallel computing with {multicore} 
#   res.prl <- runMCMC.multiChain(Y=Weed[,3], L=0, loc=Weed[,1:2], 
#                      X=NULL, MCMCinput=input.Weed, n.chn=4, n.cores=4 )

# parallel computing with {snowfall}
  res2.prl <- runMCMC.sf(Y=Weed[,3], L=0, loc=Weed[,1:2], 
                     X=NULL, MCMCinput=input.Weed, n.chn=4, n.cores=4, 
                     cluster.type="SOCK" )

#################################
### Burn-in, Thinning, and Mixing
#################################
# For the output of runMCMC()
res.m <- cutChain(res, chain.ind=1:4, burnin=500, thinning=1)

# For the output of runMCMC.multiChain() and runMCMC.sf()
res.m.prl <- lapply(res2.prl, cutChain, chain.ind=1:4, burnin=500,
                    thinning=5)
res.m <- mixChain(res.m.prl)

#################################
### Examine Posterior Samples
#################################
## Use functions in {coda} package
# convert to mcmc class
chn1.mcmc <- mcmc(cbind(beta=res.m$m, sigma=res.m$s, phi=res.m$a))
# basic information and summarized statistics
summary(chn1.mcmc)
# trace and density plots
plot(chn1.mcmc, auto.layout = TRUE)
# cross-correlation plot
crosscorr.plot(chn1.mcmc)
# auto-correlation plot
autocorr.plot(chn1.mcmc)
# effective sample size adjusted for autocorrelation
effectiveSize(chn1.mcmc)
# Gewekea??s convergence diagnostic
geweke.diag(chn1.mcmc, frac1=0.1, frac2=0.5)
# Geweke-Brooks plot
geweke.plot(chn1.mcmc, frac1=0.1, frac2=0.5)
# Heidelberger and Welcha??s convergence diagnostic
heidel.diag(chn1.mcmc, eps=0.1, pvalue=0.05)

## Plot auto-correlation curves for latent variables
plotACF(res.m$S.posterior)

## Estimate the mode of posterior samples
phi.est <- findMode(res.m$a.posterior)


######################
### Prediction
######################
## Choose unsampled locations
locp <- cbind(Weed[1:50,1]+30, Weed[1:50,2])
plot(Weed[,1:2])
points(locp, col=2, pch=2)
Lp <- rep(1, nrow(locp)) 

## Perform prediction
Ypred <- predY(res.m, loc=Weed[,1:2], locp, X=NULL, Xp=NULL, Lp=Lp, k=1, 
               rho.family="rhoPowerExp", Y.family="Poisson"
               #, parallel = "snowfall", n.cores = 4, cluster.type = "SOCK"
               ## uncomment the above line for parallel prediction
               )
Ypred.avg <- rowMeans(Ypred$Y)
plotData(Weed[,3], Weed[,1:2], Ypred.avg, locp, xlab="Eastings", ylab="Northings")


################################################
################################################
### CHAPTER 4: Model Checking
################################################
################################################

## Use Weed data and Poisson log-normal spatial model
## as an example for model checking
data(Weed)
input <- MCMCinput( run=5000, run.S=1, rho.family="rhoPowerExp", 
                    Y.family = "Poisson", ifkappa=0,
                    scales=c(0.3, 3.5, 0.9, 0.6, 0.5), 
                    phi.bound=c(0.005, 1), 
                    initials=list(c(1), 1, 0.1, 1) )
loc <- unifLoc(Weed[,1:2])
res <- runMCMC(Y=Weed[,3], L=0, loc=loc, MCMCinput=input)
res.m <- cutChain(res, chain.ind=1:4, burnin=500, thinning=5)
chn1 <- cbind(beta=res.m$m, sigma=res.m$s, phi=res.m$a)
chn1.mcmc <- mcmc(chn1)
plot(chn1.mcmc, auto.layout = TRUE)

###########################
### Bayesian Model Checking
###########################

## Define Diagnostic Statistic
# the average as diagnostic statistic
funcT <- function(Y){ mean(Y) }
# the Pearson residual type of diagnostic statistic
funcT <- function(Y){ sum((Y - mean(Y))^2/var(Y)) }

## Simulate Reference Data Sets
# Estimate parameters from posterior samples
Y.rep <- repYeb(N.sim=2000, loc, L=rep(1,nrow(loc)), res.m=res.m, est="mode")
# Pre-determined parameters
Y.rep <- repYeb(N.sim=2000, loc, L=rep(1,nrow(loc)), beta = 5, sigma = 1,
               phi = 0.1, k = 1)

## Compare Diagnostic Statistics
pRPS(T.obs = 2, T.rep = rnorm(1000))
plot_pRPS(2, rnorm(1000), nm="t")

## Bayesian model checking function
BMCT(Y.obs=Weed[,3], Y.rep, funcT, ifplot = FALSE)


#################################
### Transformed Residual Checking
#################################

## Simulate reference data sets
Y.rep <- repYeb(N.sim=5000, loc, L=rep(1,nrow(loc)), res.m=res.m, est="mode")

## Approximate transformed residuals for the observed data
etran <- tranR(Weed[,3], Y.rep, discrete = FALSE)

## Plot Transformed Residuals
plot_etran(etran, fig = 1:4)

## Calculate Hellinger Distance
d.obs <- e2dist(etran)

## Build Baseline Distribution
# simulate samples of distance
d.base <- baseline.dist(n = 100, iter = 1000) 
# baseline.parallel() (parallel version). 
# CAUSION: baseline.parallel IS REMOVED because of {multicore} issue.
# or load integrated samples
data(Dbase_n100N5000)
str(d.base)
# visualize the baseline distribution
plot_baseline(d.base[,1], colnames(d.base)[1])

## Compare the distance calculated from the observed data 
## and the corresponding baseline distribution of distance
pOne(d.obs, d.base)
plot_pRPS(d.obs[2], d.base[,2])


#################################
### END END END
#################################



