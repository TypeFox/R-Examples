"genesimulator" <-
function(nreps=3, nps=100, shape=4, scale=100){
mu <- rgamma(n=nps, shape=shape, scale=scale)
muwithreps <- rep(mu, rep(nreps, nps))
repix <- rep(1:nreps, nps)
psix <- rep(1:nps, rep(nreps,nps))
m <- cbind(muwithreps, repix, psix)
dimnames(m) <- list(NULL, c("mu", "repix", "psix"))
m

}

