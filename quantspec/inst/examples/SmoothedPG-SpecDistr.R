################################################################################
## This script illustrates how to estimate integrated quantile spectral densities

## Simulate a time series Y1,...,Y128 from the QAR(1) process discussed in
## Dette et. al (2015).
set.seed(2581)
Y <- ts1(128)

## For a defined set of quantile levels ... 
levels <- c(0.25,0.5,0.75)

## ... and a weight (of Type A), defined using the Epanechnikov kernel ...
wgt <- specDistrWeight()

## ... compute a smoothed quantile periodogram (based on the clipped time series).
## Repeat the estimation 100 times, using the moving blocks bootstrap with
## block length l=32.
sPG.cl <- smoothedPG(Y, levels.1 = levels, type="clipped", weight = wgt,
    type.boot = "mbb", B=100, l=32)

## Create a (model) spectral density kernel for he QAR(1) model for display
## in the next plot.
csd <- quantileSD(N=2^8, seed.init = 2581, type = "copula",
    ts = ts1, levels.1=levels, R = 100)
icsd <- integrQuantileSD(csd)

plot(sPG.cl, ptw.CIs = 0.1, qsd = icsd, type.CIs = "boot.full")

