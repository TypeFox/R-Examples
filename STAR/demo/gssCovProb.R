## Use a multi-core CPU and snow to perform a Monte Carlo
## investigation of the coverage probability
## of the confidence interval of gsspsth0 when
## the nominal CI is 95%

## Start by defining a function to initialize
## rstream on the different nodes
clusterSetupRSTREAM <- function(cl, ## a snow cluster
                                seed=rep(20061001,6),
                                ...) {

  nbWorkers <- length(cl)
  library(rstream)
  clusterEvalQ(cl,library(rstream))
  clusterCall(cl,
              function() {
                cmd <- parse(text=".s <- new(\"rstream.mrg32k3a\")")
                eval(cmd,env=globalenv())
              }
              )
  
  seed.txt <- paste("c(",paste(seed,collapse=","),")",sep="")
  lecuyerList <<- vector(mode="list",length=nbWorkers)
  for (wIdx in 1:nbWorkers) {
    if (wIdx == 1) cmd <- parse(text=paste("lecuyerList[[1]]",
                                  " <- new(\"rstream.mrg32k3a\",seed=",
                                  seed.txt,")",sep=""))
    else cmd <- parse(text=paste("lecuyerList[[",
                        wIdx,"]] <- new(\"rstream.mrg32k3a\")",sep=""))
    eval(cmd,env=globalenv())
    cmd <- parse(text=paste("rstream.packed(lecuyerList[[",
                   wIdx,"]]) <- TRUE",sep=""))
    eval(cmd,env=globalenv())
  }

  setLecuyer <- function(packedlecuyer) {
    assign("lecuyer",packedlecuyer,env=globalenv())
    cmd <- parse(text="rstream.packed(lecuyer)<-FALSE")
    eval(cmd,env=globalenv())
  }

  clusterApply(cl,lecuyerList,setLecuyer)
  clusterEvalQ(cl,rstream.RNG(lecuyer))

}

## load snow
library(snow)
## create a socket cluster on localhost using 2 slaves
nbSlaves <- 2
cl <- makeCluster(rep("localhost",nbSlaves),type="SOCK")
## initialize a parallel rng with independent streams on each slave
clusterSetupRSTREAM(cl)
## write a simmple function performing simulation from a
## gsspsth0 object and refiting it.
## the function returns a list of lists.
## compponent data: the simulated data
## component fit: the fitted gsspsth0 objets
quickMC2 <- function(nb, ## the number of replicates to simulate
                     gf ## the original gsspsth0 object
                     ) {
  newData <- simulate(gf,nb)
  fit <- lapply(newData, function(l) gsspsth0(l,binSize=0.025))
  list(data=newData,fit=fit)
}
## load STAR on master node
library(STAR)
## load STAR on each slave node
clusterEvalQ(cl,library(STAR))
## get data
data(e070528citronellal)
## fit the response of the first neuron
n1CitrGSSPSTH0_25 <- gsspsth0(e070528citronellal[[1]],binSize=0.025)
## simulate and fit data
## This takes 20 minutes on a "small" laptop (2x Intel Core2 2GHz, with 1 GB RAM)
nbRep <- 100
system.time(simA <- clusterApply(cl,rep(nbRep/nbSlaves,nbSlaves),quickMC2,n1CitrGSSPSTH0_25))
## Find out the number of points of the true intantaneous firing
## frequency which are out of the confidence band with a nominal
## coverage probability of 0.95
quickStats <- lapply(simA,
                     function(Bl)
                     sapply(Bl$fit,
                            function(gf) sum(n1CitrGSSPSTH0_25$freq<gf$ciLow|gf$ciUp<n1CitrGSSPSTH0_25$freq)
                            )
                     )
quickStats <- unlist(quickStats)
## make a table
table(quickStats)
## get the stats
summary(quickStats)
## compare with nominal value
length(n1CitrGSSPSTH0_25$freq)*0.05
## findout the worst case
worstIdx <- which.max(quickStats)
## plot confidence bands for worst case togther with true
## curve
idxA <- (worstIdx %/% length(simA[[1]]$fit)) + 1
idxB <- worstIdx %% length(simA[[1]]$fit)
plot(simA[[idxA]]$fit[[idxB]],colCI="grey50")
lines(n1CitrGSSPSTH0_25$mids,n1CitrGSSPSTH0_25$freq,col=2)
