## ----generate data, echo=TRUE--------------------------------------------
times <- seq(0, 46, 2)

amplitudes <- seq(0, 5, length.out=10)
means <- 1:10
data <- matrix(rnorm(24*10), nrow=10)
data <- sweep(data, 1, means, "+")

rhythms <- matrix(rep(cos(times*pi/12), 10), nrow=10, 
                  byrow=T)
rhythms <- sweep(rhythms, 1, amplitudes, "*")

rhythm.data <- data + rhythms

## ----hreg, echo=TRUE-----------------------------------------------------
library(HarmonicRegression)

hreg <- harmonic.regression(t(rhythm.data), times)

plot(amplitudes/means, hreg$pars$amp)
abline(0, 1)


