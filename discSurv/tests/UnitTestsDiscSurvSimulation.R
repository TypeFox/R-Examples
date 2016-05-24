library(discSurv)
library(Matrix)
library(matrixcalc)

#######################
# tauToPearson

# Check equal cases
stopifnot(tauToPearson (Tau=1) == 1)
stopifnot(tauToPearson (Tau=0) == 0)
stopifnot(tauToPearson (Tau=-1) == -1)

# Check if Function is symmetric
stopifnot(all.equal(-sapply(seq(-1, 0, length.out=100), tauToPearson), sapply(seq(1, 0, length.out=100), tauToPearson)))

########################
# Design of Simulation

# Monte carlo samples for each design = 10
# SampleSize = 100 in each replicate

# Responses: 
# R1 i.n.d. ~ F (df1=3, df2=5, ncp=exp(\eta))
# R2 i.n.d. ~ F (df1=5, df2=3, ncp=exp(-\eta))
# R3 i.n.d. ~ F (df1=3, df2=3, ncp=exp(\eta))

# Random Censoring
# Independent of survival times C_i i.i.d ~ Gamma (shape=1, scale=2)

# Correlation structure of survival times:
# Specify with kendalls tau -> spearmans rho
# Kendalls tau matrix:
#     R1   R2   R3
# R1  1   0.3  0.4
# R2  0.3   1  0.5
# R3  0.4 0.5    1

# Covariates: 
# V1: Binary variable ~ Bin (n=4, \pi=0.25) -> E (V1) = 1
# V2: Continuous positive variable ~ Gamma (shape=1, scale=1) -> E (V2) = 1
# V3: Continuous variable ~ Normal (\mu=0, \sigma^2=1) -> E (V5) = 0

# True linear predictor:
# \eta = \boldsymbol {X} %*% \beta
# \beta_1 = c(-0.5, 1, 0.5)
# \beta_2 = c(1, -1, 1)
# \beta_3 = c(-0.5, -0.5, 2)

# Correlation Structure in Covariates:
# Specify with kendalls tau -> pearsons rho
#         V1    V2     V3
# V1       1 -0.25      0
# V2   -0.25     1   0.25
# V3       0  0.25      1

# Design Correlation Matrix of covariates
DesignCovariateCor <- diag(3)
DesignCovariateCor [lower.tri(DesignCovariateCor)] <- c(-0.25, 0, 0.25)
DesignCovariateCor [upper.tri(DesignCovariateCor)] <- c(-0.25, 0, 0.25)
# Check if symmetric
sum(DesignCovariateCor-t(DesignCovariateCor))==0 # TRUE -> ok
# Check if positive definite
is.positive.definite (DesignCovariateCor) # TRUE -> ok
# Check if transformed pearson matrix is positive, definite
is.positive.definite (apply(DesignCovariateCor, c(1,2), tauToPearson)) # TRUE -> ok

# Design Correlation Matrix positive definite after transformation
DesignResponseCor <- diag(3)
DesignResponseCor [lower.tri(DesignResponseCor)] <- c(0.3, 0.4, 0.5)
DesignResponseCor [upper.tri(DesignResponseCor)] <- c(0.3, 0.4, 0.5) 

# Check if symmetric
sum(DesignResponseCor-t(DesignResponseCor))==0
# Check if positive definite
is.positive.definite (DesignResponseCor)
# Check if transformed pearson matrix is positive, definite
is.positive.definite (apply(DesignResponseCor, c(1,2), tauToPearson))

# Simulate raw data
DesignSampleSize <- 100
MonteCarloSamples <- 10
SimOutput <- vector("list", MonteCarloSamples)
for(j in 1:MonteCarloSamples) {
  SimOutput [[j]] <- simCompRisk (responseCorr=DesignResponseCor, covariateCorr=DesignCovariateCor, covariateSeed=NULL, sampleSize=100, 
                                  covariateQuantFunc=c("qbinom", "qgamma", "qnorm"), covariateArgs=list(c(size=4, prob=0.25), c(shape=1, scale=1), c(mean=0, sd=1)), 
                                  intercept=c(FALSE, FALSE, FALSE), trueCoef=list(c(-0.5, 1, 0.5), c(1, -1, 1), c(-0.5, -0.5, 2)), 
                                  responseSeed=NULL, responseQuantFunc=c("qf", "qf", "qf"), responseFixArgs=list(c(df1=3, df2=5), c(df1=5, df2=3), c(df1=3, df2=3)), responseLinPredArgs=list(list(ncp=function (x) {exp(x)}), list(ncp=function (x) {exp(-x)}), list(ncp=function (x) {exp(x)})), 
                                  censorRN="rgamma", censorArgs=c(shape=1, scale=2), censorSeed=NULL)
}
SimOutput [[1]]
SimOutput [[5]]
SimOutput [[10]]

OverviewData <- data.frame(Min=rep(NA, MonteCarloSamples), Max=rep(NA, MonteCarloSamples), Censrate=rep(NA, MonteCarloSamples))
for(j in 1:MonteCarloSamples) {
  
  # Calculate Range of observed time
  OverviewData [j, 1:2] <- range (SimOutput [[j]]$Data [, "Time"])
  
  # Calculate censoring rate
  OverviewData [j, "Censrate"] <- mean(SimOutput [[j]]$Data [, "Censor"])
  
}
OverviewData

# Additional checks
for(i in 1:MonteCarloSamples) {
  stopifnot(all(sapply(1:4, function (x) all(SimOutput [[i]]$Data [, 2:5] [, x] == 1 | SimOutput [[i]]$Data [, 2:5] [, x] == 0))))
  stopifnot(all(SimOutput [[i]]$Data [, 1]>=0))
}
stopifnot(all(OverviewData [, 3] >= 0))
