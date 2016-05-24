################################################################################
## This script illustrates how to work with QuantilePG objects

## Simulate a time series Y1,...,Y128 from the QAR(1) process discussed in
## Dette et. al (2015).
Y <- ts1(64)

## For a defined set of quantile levels
levels <- c(0.25,0.5,0.75)

## the various quantile periodograms can be calculated calling quantilePG:

## For a copula periodogram as in Dette et. al (2015) the option 'type="qr"'
## has to be used:
system.time(
    qPG.qr <- quantilePG(Y, levels.1 = levels, type="qr"))

## For the CR-periodogram as in Kley et. al (2016) the option 'type="clipped"'
## has to be used. If bootstrap estimates are to be used the parameters
## type.boot, B and l need to be specified.
system.time(
    qPG.cl <- quantilePG(Y, levels.1 = levels, type="clipped",
        type.boot="mbb", B=250, l=2^5))

## The two previous calls also illustrate that computation of the CR-periodogram
## is much more efficient than the quantile-regression based copula periodogram.

## Either periodogram can be plotted using the plot command
plot(qPG.cl)
plot(qPG.qr)

## Because the indicators are not centered it is often desired to exclude the
## frequency 0; further more the frequencies (pi,2pi) are not wanted to be
## included in the plot, because f(w) = Conj(f(2 pi - w)).
## Using the plot command it is possible to select frequencies and levels for
## the diagram:
plot(qPG.cl, frequencies=2*pi*(1:32)/64, levels=c(0.25))

## We can also plot the same plot together with a (simulated) quantile spectral
## density kernel
csd <- quantileSD(N=2^8, seed.init = 2581, type = "copula",
    ts = ts1, levels.1=c(0.25), R = 100)
plot(qPG.cl, qsd = csd, frequencies=2*pi*(1:32)/64, levels=c(0.25))

## Calling the getValues method allows for comparing the two quantile
## periodograms; here in a diagram:
freq <- 2*pi*(1:31)/32
V.cl <- getValues(qPG.cl, frequencies = freq, levels.1=c(0.25))
V.qr <- getValues(qPG.qr, frequencies = freq, levels.1=c(0.25))
plot(x = freq/(2*pi), Re(V.cl[,1,1,1]), type="l",
        ylab="real part -- quantile PGs", xlab=expression(omega/2*pi))
lines(x = freq/(2*pi), Re(V.qr[,1,1,1]), col="red")

## Now plot the imaginary parts of the quantile spectra for tau1 = 0.25
## and tau2 = 0.5
freq <- 2*pi*(1:31)/32
V.cl <- getValues(qPG.cl, frequencies = freq, levels.1=c(0.25, 0.5))
V.qr <- getValues(qPG.qr, frequencies = freq, levels.1=c(0.25, 0.5))
plot(x = freq/(2*pi), Im(V.cl[,1,2,1]), type="l",
    ylab="imaginary part -- quantile PGs", xlab=expression(omega/2*pi))
lines(x = freq/(2*pi), Im(V.qr[,1,2,1]), col="red")
