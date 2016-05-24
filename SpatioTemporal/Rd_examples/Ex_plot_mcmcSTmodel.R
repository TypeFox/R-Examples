##load MCMC results instead
data(MCMC.mesa.model)

##plot the log-likelihood
plot(MCMC.mesa.model, ylab="", xlab="", type="l")

##and MCMC tracks for four of the parameters
par(mfrow=c(4,1),mar=c(2,2,2.5,.5))
for(i in c(4,9,13,15)){
  plot(MCMC.mesa.model, i, ylab="", xlab="", type="l")
}

##or by name
par(mfrow=c(1,1),mar=c(2,2,2.5,.5))
plot(MCMC.mesa.model, "nu.log.range.exp", ylab="", xlab="", type="l",
     main="all range estimates", ylim=c(-14,10))
##all ranges in one plot
plot(MCMC.mesa.model, "log.range.const.exp", add=TRUE, col=2)
plot(MCMC.mesa.model, "log.range.V1.exp", add=TRUE, col=3)
plot(MCMC.mesa.model, "log.range.V2.exp", add=TRUE, col=4)
