##load estimation results
data(est.mesa.model)
##and MCMC results instead
data(MCMC.mesa.model)

##compute density estimates for the results, and use the Gaussian approximation
##based on Fischer information as reference.
dens <- density(MCMC.mesa.model, estSTmodel=est.mesa.model)

##all the estimated densities
str(dens,1)

##or results for one paramter
dens[[1]]

##plot density functions
plot(dens)
##for a different paramter, along with Gaussian approx
plot(dens, 3, norm.col="red")

##all covariance parameters
par(mfrow=c(3,3),mar=c(4,4,2.5,.5))
for(i in 9:17){
  plot(dens, i, norm.col="red")
}
