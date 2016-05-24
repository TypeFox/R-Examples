## This script can be used to create and store a QuantileSD object

## Parameters for the simulation:
R <- 50                      # number of independent repetitions;
                             # R should be much larger than this in practice!
N <- 2^8                     # number of Fourier frequencies in [0,2pi)
ts <- ts1                    # time series model
levels <- seq(0.1,0.9,0.1)   # quantile levels
type <- "copula"             # copula, not Laplace, spectral density kernel
seed.init <- 2581            # seed for the pseudo random numbers

## Simulation takes place once the constructor is invoked
qsd <- quantileSD(N=N, seed.init = 2581, type = type,
    ts = ts, levels.1=levels, R = R)

## The simulated copula spectral density kernel can be called via
V1 <- getValues(qsd)

## It is also possible to fetch the result for only a few levels
levels.few <- c(0.2,0.5,0.7)
V2 <- getValues(qsd, levels.1=levels.few, levels.2=levels.few)

## If desired additional repetitions can be performed to yield a more precise
## simulation result by calling; here the number of independent runs is doubled.
qsd <- increasePrecision(qsd,R)

## Often the result will be stored for later usage.  
save(qsd, file="QAR1.rdata")

## Take a brief look at the result of the simulation
plot(qsd, levels=levels.few)

## When plotting more than only few levels it may be a good idea to plot to
## another device; e. g., a pdf-file
K <- length(levels)
pdf("QAR1.pdf", width=2*K, height=2*K)
  plot(qsd)
dev.off()

## Now we analyse the multivariate process (eps_t, eps_{t-1}) from the
## introduction of Barunik&Kley (2015). It can be defined as
ts_mult <- function(n) {
  eps <- rnorm(n+1)
  return(matrix(c(eps[2:(n+1)], eps[1:n]), ncol=2))
}

## now we determine the quantile cross-spectral densities
qsd <- quantileSD(N=N, seed.init = 2581, type = type,
    ts = ts_mult, levels.1=levels, R = R)

## from which we can for example extract the quantile coherency
Coh <- getCoherency(qsd, freq = 2*pi*(0:64)/128)

## We now plot the real part of the quantile coherency for j1 = 1, j2 = 2,
## tau1 = 0.3 and tau2 = 0.6
plot(x = 2*pi*(0:64)/128, Re(Coh[,1,3,2,6]), type="l")
