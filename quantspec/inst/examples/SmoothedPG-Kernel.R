################################################################################
## This script illustrates how to work with QuantilePG objects

## Simulate a time series Y1,...,Y128 from the QAR(1) process discussed in
## Dette et. al (2015).
Y <- ts1(128)

## For a defined set of quantile levels ... 
levels <- c(0.25,0.5,0.75)
FF <- 2*pi*(0:64)/128
K <- length(levels)

## ... and a weight (of Type A), defined using the Epanechnikov kernel ...
wgt <- kernelWeight(W=W1, N=128)

## ... compute a smoothed quantile periodogram (based on the clipped time series).
## Repeat the estimation 100 times, using the moving blocks bootstrap with
## block length l=32.
sPG.cl <- smoothedPG(Y, levels.1 = levels, type="clipped", weight = wgt,
    type.boot = "mbb", B=100, l=32)

## There is a getValues command that works analogously to the one for
## objects of type QuantilePG 
V <- getValues(sPG.cl, frequencies = 2*pi*(1:31)/32,
    levels.1=c(0.25), levels.2=c(0.5,0.75))

## Now create three plots with the different types of confidence intervals.
plot(sPG.cl, type.CIs = "naive.sd")
plot(sPG.cl, type.CIs = "boot.sd")
plot(sPG.cl, type.CIs = "boot.full")


## Create a (model) spectral density kernel for he QAR(1) model for display
## in the next plot.
csd <- quantileSD(N=2^8, seed.init = 2581, type = "copula",
    ts = ts1, levels.1=levels, R = 100)

## Now show the same plot, together with the (simulated) copula spectral
## density kernel and the periodogram that was used for smoothing;
## for later usage we create the plot in a pdf-file.
pdf("plot.pdf", width=3*K, height=3*K)
  plot(sPG.cl, plotPG = TRUE, qsd = csd,
    frequencies = FF[FF > 0], type.CIs = "naive.sd", type.scaling="individual")
dev.off()

## Note that various combinations are possible:
plot(sPG.cl, plotPG = TRUE, ptw.CIs = 0)
plot(sPG.cl, plotPG = TRUE, frequencies = FF[FF > 0], ptw.CIs = 0)
plot(sPG.cl, plotPG = TRUE, frequencies = FF[FF > 0], ptw.CIs = 0.1)
plot(sPG.cl, plotPG = FALSE, ptw.CIs = 0.1, type.CIs = "boot.sd")
plot(sPG.cl, plotPG = FALSE, qsd = csd, ptw.CIs = 0.1, type.CIs = "boot.sd")
plot(sPG.cl, plotPG = TRUE, frequencies = FF[FF > 0], ptw.CIs = 0, levels=c(0.25))

## and many more...