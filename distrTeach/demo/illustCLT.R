## illustrateCLT produces plots for the subseqent laws of T_n
## syntax: illustrateCLT(Distr, len, sleep = 0) where 
##         Distr ^= distribution of X_i, 
##           len ^= # summands at which to stop
##         sleep ^= duration of the pause between the different plots

require(distrTeach)
options("newDevice"=TRUE)
options("device.ask.default"=FALSE)
# some examples
# distroptions("DefaultNrFFTGridPointsExponent" = 13)
illustrateCLT(Distr = Unif(), len = 20, sleep = 0.5)
illustrateCLT(Distr = Exp(), len = 20, sleep = 0.5)

N <- Norm(mean = 2, sd = 1.3)
P <- Pois(lambda = 1.2)
Z <- 2 * N + 3 + P # exact transformation
illustrateCLT(Distr = sin(abs(Z)), len = 20, sleep = 0.5) #something weird

#illustrateCLT(Distr = Chisq(), len = 20)
#illustrateCLT(Distr = Td(df = 5), len = 20)
#illustrateCLT(Distr = Beta(), len = 20)
#distroptions("DefaultNrFFTGridPointsExponent" =  14)
#illustrateCLT(Distr = Lnorm(), len = 20)

distroptions("DefaultNrFFTGridPointsExponent" = 12)
illustrateCLT(Distr = Pois(lambda = 2), len = 20, sleep = 0.5)
illustrateCLT(Distr = Binom(size = 5), len = 20, sleep = 0.5)
illustrateCLT(Distr = Nbinom(size = 5), len = 20, sleep = 0.5)
