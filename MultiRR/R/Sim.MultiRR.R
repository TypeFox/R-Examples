Sim.MultiRR <-
function(n.ind, SeriesPerInd, ObsPerLevel, EnvGradient, PopInt, PopSlope, VCVInd, VCVSeries, ResVar, n.sim, unbalanced=FALSE, prop.ind, complete.observations=TRUE, n.obs) {
SimData <- list()
ss <- expand.grid(n.ind, SeriesPerInd)
colnames(ss) <- c("Individuals", "SeriesPerIndividual")
DataFrames <- array(list(),c(n.sim, nrow(ss)))
IndividualsD <- array(NA,c(n.sim, nrow(ss)))
SeriesD <- array(NA,c(n.sim, nrow(ss)))

for(j in 1:nrow(ss)){
n.ind <-ss[j,1]
SeriesPerInd <-ss[j,2]

for(i in 1:n.sim){    
if (unbalanced==TRUE) {
SeriesPerInd2 <- sample(1:SeriesPerInd, size=n.ind, replace=TRUE, prob=prop.ind)
IDSeries <- rep(1:n.ind, SeriesPerInd2)
TotalSeries <- length(IDSeries)
repsID <- SeriesPerInd2*ObsPerLevel
} else {
SeriesPerInd2 <- SeriesPerInd
TotalSeries <- n.ind*SeriesPerInd2
repsID <- rep(SeriesPerInd2*ObsPerLevel, n.ind)
}

obs <- length(EnvGradient)*TotalSeries
BetaInd <- mvrnorm(n.ind, c(0,0), VCVInd)
IntInd <- BetaInd[,1]
SlopeInd <- BetaInd[,2]

BetaSeries <- mvrnorm(TotalSeries, c(0,0), VCVSeries)
IntSeries <- BetaSeries[,1]
SlopeSeries <- BetaSeries[,2]

DataFrame.1 <- data.frame(PopInt=rep(PopInt, obs), PopSlope=rep(PopSlope, obs), Ind=as.factor(rep(seq(1:n.ind), repsID)), IntInd=rep(IntInd, repsID), SlopeInd=rep(SlopeInd, repsID), Series=as.factor(rep(seq(1:TotalSeries), each=length(EnvGradient))), IntSeries=rep(IntSeries, each=length(EnvGradient)), SlopeSeries=rep(SlopeSeries, each=length(EnvGradient)), x=EnvGradient)

DataFrame.2 <- DataFrame.1[rep(seq_len(nrow(DataFrame.1)), ObsPerLevel), ]


if (complete.observations==TRUE) {n.obs=nrow(DataFrame.2)
              } else {
                  n.obs=n.obs
              }

DataFrame <-  DataFrame.2[sample(1:nrow(DataFrame.2), n.obs),]


DataFrame$y=(DataFrame$PopInt + DataFrame$IntInd + DataFrame$IntSeries) + ((DataFrame$PopSlope + DataFrame$SlopeInd + DataFrame$SlopeSeries)*DataFrame$x) + rnorm(nrow(DataFrame), 0, sqrt(ResVar))

DataFrames[[i,j]] <-DataFrame 
IndividualsD[[i,j]] <- length(levels(DataFrame$Ind))
SeriesD[[i,j]] <- length(levels(DataFrame$Series))

}
}
SimData$DataFrames <- DataFrames
SimData$VCVInd <- VCVInd
SimData$VCVSeries <- VCVSeries
SimData$n.sim <- n.sim
SimData$Individuals <- IndividualsD
SimData$SeriesD <-SeriesD 
SimData$SimInt <- PopInt
SimData$SimSlope <- PopSlope
SimData$Residuals <- ResVar
SimData
}
