## illustrateLLN produces plots for the subseqent laws of sum X_i/n
## syntax:   illustrateLLN <- function(Distr, n, m, step, sleep, ...)
##         Distr ^= distribution of X_i, 
##                  default: Norm()   
##             n ^= vector of sample sizes to be considered 
##                  default: c(1,3,5,10,25,50,100,500,1000,10000)
##             m ^= # replicates of sum X_i/n to be plotted
##                  default: 50
##          step ^= # how many (new) replicates to be plotted in one step
##                  default: 1
##         sleep ^= duration of the pause between the different plots
##                  default: 0
##           ... ^= further arguments to plot

require(distrTeach)
options("newDevice"=TRUE)
options("device.ask.default"=FALSE)
# some examples
# distroptions("DefaultNrFFTGridPointsExponent" = 13)
illustrateLLN(Distr = Norm(0,3), sleep = 0.1)
illustrateLLN(Distr = Unif(), sleep = 0.1)
illustrateLLN(Distr = Exp(), sleep = 0.1)

N <- Norm(mean = 2, sd = 1.3)
P <- Pois(lambda = 1.2)
Z <- 2 * N + 3 + P # exact transformation
illustrateLLN(Distr = sin(abs(Z)), sleep = 0.1) #something weird

illustrateLLN(Distr = Pois(lambda = 2), sleep = 0.1)
illustrateLLN(Distr = Binom(size = 5), sleep = 0.1)
illustrateLLN(Distr = Nbinom(size = 5), sleep = 0.1)
