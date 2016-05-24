library(rjags)
library(covafillr)

loadJAGSModule()

fun <- function(x) x ^ 2
n <- 100
x <- runif(n,-2,2)
y <- rnorm(n,fun(x),0.5)
obsC <- seq(-3,3,len=1000)
obs <- fun(obsC) + rnorm(length(obsC),0,0.1)


plot(obsC,obs)


jags <- jags.model('covafill.jags',
                   data = list(N = n,
                               x = matrix(x,ncol=1),
                               y = y,
                               obsC = matrix(obsC,ncol=1),
                               obs = obs,
                               h = c(1)),
                   n.chains = 1,
                   n.adapt = 100)

samp <- jags.samples(jags,c('sigma','cf'),n.iter = 10000, thin = 10)

plot(x,samp$cf[,1,1])
hist(samp$sigma)


## 2D

fun <- function(x) 5 * x[1]^2
n <- 10000
x <- cbind(runif(n,-2,2), runif(n,5,7))
y <- rnorm(n,apply(x,1,fun),0.5)
obsC <- expand.grid(seq(-3,3,len=100),seq(4,8,len=100))
obs <- apply(obsC,1,fun) + rnorm(dim(obsC)[1],0,0.1)


lattice::levelplot(obs ~ obsC[,1] + obsC[,2])


load.module('covafillr',system.file('libs',package='covafillr'))

jags <- jags.model('covafill.jags',
                   data = list(N = n,
                               x = x,
                               y = y,
                               obsC = obsC,
                               obs = obs,
                               h = c(0.5,0.5)),
                   n.chains = 1,
                   n.adapt = 100)

samp <- jags.samples(jags,c('sigma','cf'),n.iter = 10000, thin = 10)

par(mfrow=c(1,2))
plot(x[,1],samp$cf[,1,1])
plot(x[,2],samp$cf[,1,1])
