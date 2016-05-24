################################################################################
# Analysis of the European Wheat Price Time Series:
# -------------------------------------------------
# In this demo the S&P 500 time series is analysed ans some features are
# explained.
################################################################################

# Take a look at the data first:
plot(wheatprices)

# Define the quantile orders.
taus <- c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)

# Determine the rank-based Laplace periodogram kernel
# for the quantile orders taus.
LPG <- quantilePG(wheatprices, levels.1 = taus, type="qr")

# Then, plot for omega in (0,pi] and various levels
f <- getFrequencies(LPG)
plot(LPG, levels=c(0.1,0.5), frequencies=f[which(f > 0 & f <= pi)])
plot(LPG, frequencies = f[which(f > 0 & f <= pi)], levels = c(0.05, 0.5, 0.95))
plot(LPG, frequencies = f[which(f > 0 & f <= pi)], levels = c(0.25, 0.5, 0.75))

# Now determine the smoothed rank-based Laplace periodogram kernel
# (from the unsmoothed one), using a weight determined by an
# Epanechnikov kernel and b=0.05
sLPG <- smoothedPG(LPG, weight=kernelWeight(W=W1, bw=0.05))

# Then, plot the estimate, with a full bootstrap pointwise confidence band
plot(sLPG, levels=c(0.25,0.5,0.75), ptw.CIs=0)

# Note the many peaks between 0.06 and 0.07, which correspond to the
# well known cycle of 15.3 years in the data.
# (the peak is at the reciprocal of 15.3!

# Determine the CR periodogram kernel
# for the quantile orders taus.
CR <- quantilePG(wheatprices, levels.1 = taus, type="clipped")

# Then, plot for omega in (0,pi] and various levels
f <- getFrequencies(CR)
plot(CR, frequencies = f[which(f > 0 & f <= pi)], levels = c(0.25, 0.5, 0.75))

# Now determine the smoothed CR periodogram kernel (from the unsmoothed one),
# using a weight determined by an Epanechnikov kernel and b=0.05
sPG <- smoothedPG(CR, weight=kernelWeight(W=W1, bw=0.05))

# plot for omega in [0,pi],
#  - for tau_1, tau_2 = 0.25, 0.5, 0.75
#  - with each subplot scaled individually to their min and max values,
#  - with pointwise (1-0.1)-confidence bands,
#  - with the unsmoothed periodogram ommitted, and
#  - for omega in [0,pi]
freq <- getFrequencies(sPG)
plot(sPG, levels=c(0.25,0.5,0.75), type.scaling="individual",
    ptw.CIs = 0.1, plotPG=F,
    frequencies=freq[which(freq >= 0 & freq <= pi)])

# Now determine the smoothed rank-based Laplace periodogram,
#  - directly from the time series,
#  - for tau_1, tau_2 in taus,
#  - also determine 250 bootstrap replicates
#    (use moving blocks and blocklength 32),
sPG <- smoothedPG(wheatprices, levels.1 = taus, weight=kernelWeight(W=W1, bw=0.05),
    type = "clipped", type.boot="mbb", B=250, l=32)

# Plot (for the values 0.25, 0.5 and 0.75),
# with homogeneous range in real and imaginary parts
plot(sPG, levels=c(0.25,0.5,0.75), type.scaling="real-imaginary")

# ... and with the pointwise confidence bands determined from the
# bootstrap replicates.
plot(sPG, levels=c(0.25,0.5,0.75), type.CIs = "boot.full")
