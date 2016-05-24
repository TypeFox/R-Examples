################################################################################
# Simulation study, analyzing the QAR(1) model from Dette et. al (2014+):
# -----------------------------------------------------------------------
# In this demo a number of R quantile periodograms and smoothed quantile
# periodograms are computed from R independent, simulated time series from the
# QAR(1) model with parameters as in the above paper.
#
# The root integrated mean squared errors are stored and can be used for
# comparison of the estimators.
#
# Finally, plots for the last simulation run are shown.
################################################################################

ts <- ts1
type <- "copula"

N <- 128
R <- 100

freq <- 2*pi*(1:16)/32
levels <- c(0.25, 0.5, 0.75)

J <- length(freq)
K <- length(levels)

# First determine the copula spectral density
csd <- quantileSD(N=2^8, seed.init = 2581, type = type,
                  ts = ts, levels.1=levels, R = 100)

# init an array for the root integrated mean squared errors
# (1: CR periodogram, 2: rank-based Laplace periodogram,
#  3: smoothed CR periodogram, 4: smoothed rank-based Laplace pg.)
sims  <- array(0, dim=c(4,R,J,K,K))

weight <- kernelWeight(W=W1, bw=0.3)

for (i in 1:R) {
  Y <- ts1(N)

  CR <- quantilePG(Y, levels.1=levels, type="clipped")
  sims[1,i,,,] <- getValues(CR, frequencies=freq)[,,,1]

  LP <- quantilePG(Y, levels.1=levels, type="qr")
  sims[2,i,,,] <- getValues(LP, frequencies=freq)[,,,1]

  sCR <- smoothedPG(CR, weight=weight)
  sims[3,i,,,] <- getValues(sCR, frequencies=freq)[,,,1]

  sLP <- smoothedPG(LP, weight=weight)
  sims[4,i,,,] <- getValues(sLP, frequencies=freq)[,,,1]
}

trueV <- getValues(csd, frequencies=freq)
SqDev <- array(apply(sims, c(1,2),
        function(x) {abs(x-trueV)^2}), dim=c(J,K,K,4,R))
rimse <- sqrt(apply(SqDev, c(2,3,4), mean))

# Inspect the rimse, but note that to have reliable approximations of
# the true rimse the simulations should run much longer.
rimse

# Finally take a look a the last simulated periodograms and see that
# time series of length 64 are not quiet long enought to yield
# reliable estimators.

f <- getFrequencies(sCR)
plot(sCR, qsd=csd, frequencies = f[f > 0])
plot(sCR, qsd=csd, plotPG=TRUE, frequencies = f[f > 0])
plot(sLP, qsd=csd, frequencies = f[f > 0],
    ptw.CIs = 0, type.scaling="real-imaginary")
plot(sLP, qsd=csd, plotPG=TRUE, frequencies = f[f > 0],
    ptw.CIs = 0, type.scaling="real-imaginary")
