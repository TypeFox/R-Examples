################################################################################
# Analysis of the S&P 500 stock index, 2004-2006:
# -----------------------------------------------
# In this demo the S&P 500 time series is analyzed and some features are
# explained.
################################################################################

require(zoo)

# Take a look at the data first:
plot(sp500)

# Ignore that the observations are not equally spaced
Y <- coredata(sp500)

# There is barely any correlation in the data
acf(Y)
spectrum(Y, method="ar", log="no", ylim=c(0, 5e-4))
spectrum(Y, method="pgram", span=75, log="no", ylim=c(0, 5e-4))

# In the series of squared returns significant autocorrelation
# can be seen. 
acf(Y^2)

# Now determine CR periodgram kernels
taus <- c(0.05,0.1,0.5,0.9,0.95)
CR <- quantilePG(Y, levels.1 = taus, type="clipped",
    type.boot="mbb", B=250, l=32)
freq <- getFrequencies(CR)
plot(CR, levels=c(0.05,0.5), frequencies=freq[which(freq > 0 & freq <= pi)])

# From the CR periodogram kernels compute a smoothed version
# use a Epanechnikov kernel and bandwidth of 0.07
sPG <- smoothedPG(CR, weight=kernelWeight(W=W1, bw=0.07))

# For the upcoming plots use only frequencies from (0,pi).
f <- freq[which(freq > 0 & freq <= pi)]

# In the extreme quantiles a peak at low frequencies is visible.
plot(sPG, levels=c(0.05,0.5,0.95), type.scaling="individual",
    frequencies=f, ptw.CIs = 0.1, plotPG=F)

# This becomes even more evident when the same scaling for
# real and imaginary parts is used;
# for a change plot confidence intervals determined from the
# bootstrapped estimator. 
plot(sPG, levels=c(0.05,0.5,0.95), type.scaling="real-imaginary",
    frequencies=f, ptw.CIs = 0.1, type.CIs="boot.full")

# We can also estimate the spectra using a lag-window type estimator,
# based on Clipped Covariances. The main difference is that smoothing
# is applied in the time domain. 
lagOp <- clippedCov(Y,levels.1 = c(0.05,0.5,0.95))
weight <- lagKernelWeight(W = WParzen,  bw = 25, K = length(Y))
lagEst <- lagEstimator(lagOp,weight = weight) 
plot(lagEst,ptw.CIs = 0.1,levels = c(0.05,0.5,0.95),type.scaling= "individual",type.CIs = "naive.sd")


