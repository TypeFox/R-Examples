### R code from vignette source 'MCMCintro.Rnw'

###################################################
### code chunk number 1: MCMCintro.Rnw:34-42
###################################################
minmaxpost <- function(theta, data){
  mu <- theta[1]
  sigma <- exp(theta[2])
  dnorm(data$min, mu, sigma, log=TRUE) +
    dnorm(data$max, mu, sigma, log=TRUE) +
    (data$n - 2) * log(pnorm(data$max, mu, sigma) -
    pnorm(data$min, mu, sigma))
}


###################################################
### code chunk number 2: MCMCintro.Rnw:51-55
###################################################
data <- list(n=10, min=52, max=84)
library(LearnBayes)
fit <- laplace(minmaxpost, c(70, 2), data)
fit


###################################################
### code chunk number 3: MCMCintro.Rnw:60-64
###################################################
mycontour(minmaxpost, c(45, 95, 1.5, 4), data,
          xlab=expression(mu), ylab=expression(paste("log ",sigma)))
mycontour(lbinorm, c(45, 95, 1.5, 4),
          list(m=fit$mode, v=fit$var), add=TRUE, col="red")


###################################################
### code chunk number 4: MCMCintro.Rnw:73-78
###################################################
mcmc.fit <-  rwmetrop(minmaxpost, 
             list(var=fit$v, scale=3), 
             c(70, 2), 
             10000, 
             data)


###################################################
### code chunk number 5: MCMCintro.Rnw:82-83
###################################################
mcmc.fit$accept


###################################################
### code chunk number 6: MCMCintro.Rnw:88-92
###################################################
mycontour(minmaxpost, c(45, 95, 1.5, 4), data,
          xlab=expression(mu), 
          ylab=expression(paste("log ",sigma)))
points(mcmc.fit$par)


###################################################
### code chunk number 7: MCMCintro.Rnw:105-110
###################################################
mu <- mcmc.fit$par[, 1]
sigma <- exp(mcmc.fit$par[, 2])
P.75 <- mu + 0.674 * sigma
plot(density(P.75), 
     main="Posterior Density of Upper Quartile")


