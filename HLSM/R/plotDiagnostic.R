
##get diagnostic plot of the chain of interest
##chains can be extracted from get' ' funcions and passed to plot.HLSM
plotDiagnostic = function(chain){
    par(mfrow = c(2,2))
    draws = length(chain)
    plot(1:draws,chain,type='l', xlab = 'Iterations', ylab = 'Parameter', main = 'Trace plot')
    ###Running Means####
    rmean = runningmeans(chain, 1)
    m=0.5
    plot(c(1:draws), rmean, type="l", xlab="", ylab="Mean", main = 'Running means')
    ##Autcorr plot##
    bozo=autocorr(mcmc(chain), lags=seq(1,draws/2,50))
    plot(seq(1,draws/2,50), bozo, type="h", ylim=c(-1,1), xlab="", ylab="", main = 'Autocorrelation')
}
