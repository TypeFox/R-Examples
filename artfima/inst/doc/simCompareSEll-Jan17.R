#Source: simCompareSE-ll-Jan17.R
#Purpose: artfima parameter estimates & observed se
#requires: arfima, ltsa, gsl
#
#Set directory and output file here
#setwd("D:/R/2015/artfima/simulations/Dec14") #home
setwd("E:/dropbox/R/2015/artfima/simulations/Jan19") #office
#Set NSIM and file output name
#NSIM <- 100 #RStudio: ; R 231 sec
NSIM <- 10^4 #R  24758 sec, 6.9 hr

#R workspace saved here with unique FILEOUT
FILEOUT <- paste0("WS-Jan17-ll",NSIM,".RData")
#
require("artfima")
require("parallel")
NCORES <- 8
#
#OneIt - modified RGN for parallel
OneIt <- function(x){
  Ns <- c(200, 500, 1000, 5000)
  out <- numeric(0)
  LAM0 <- x$LAM
  D0 <- x$D
  for (jn in 1:length(Ns)) {
    n <- Ns[jn]
    z <- artSim(n=n, lambda=LAM0, d=D0)
    ansW <- artfima(z, likAlg="Whittle")
    outj <- c(n, LAM0, D0, ansW$lambdaHat, ansW$dHat, ansW$se, ansW$convergence,
              as.numeric(ansW$onBoundary)) #has length 9 since se has 2 components
    out <- c(out, outj)
  }
  out
}
#set define parameters for the simulation
Ns <- c(200,500,1000,5000)#copy defined in OneIt. Needed for parallel.
LAMs <- c(0.01, 0.02, 0.03)
Ds <- c(0.4, 0.8, 1.2)
NumPars <- length(LAMs)*length(Ds)
#set w matrix
w <- vector(mode="list", length=NumPars*NSIM)
it <- 0
for (k in 1:length(LAMs)) {
  for (i in 1:length(Ds)) {
    for (jsim in 1:NSIM) {
      it <- it+1
      w[[it]] <- list(LAM=LAMs[k], D=Ds[i])
    }
  }
}
#using parLapply
date()
startTime <- proc.time()[3]
cl <- makeCluster(spec=NCORES, type="PSOCK")
clusterSetRNGStream(cl, iseed=775123)
#initialized seed to a fixed default
#Export variables
clusterExport(cl, list("OneIt"))
#Export library
clusterEvalQ(cl, require("artfima"))
#parallel lapply
OUT<-parLapply(cl, w, fun=OneIt)
stopCluster(cl) #stop cluster
date()
totalTime <- proc.time()[3]-startTime
totalTime
save.image(file=FILEOUT)

