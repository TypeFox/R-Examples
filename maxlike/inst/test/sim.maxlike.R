


# Simulation example

sim <- function() {
    x1 <- sort(rnorm(100))
    x1 <- raster(outer(x1, x1), xmn=0, xmx=100, ymn=0, ymx=100)
    x2 <- raster(matrix(runif(1e4), 100, 100), 0, 100, 0, 100)

    logit.psi <- -1 + 1*x1 + 0.5*x2
    psi <- exp(logit.psi)/(1+exp(logit.psi))

    r <- stack(x1, x2)
    r@layernames <- c("x1", "x2")

    pa <- matrix(NA, 100, 100)
    pa[] <- rbinom(1e4, 1, as.matrix(psi))
    pa <- raster(pa, 0, 100, 0, 100)
    xy <- xyFromCell(pa, sample(Which(pa==1, cells=TRUE), 1000))

    return(list(r=r, xy=xy))
}

str(sim(), 1)

nsim <- 100
simout <- matrix(NA, nsim, 3)
colnames(simout) <- c('beta0', 'beta1', 'beta2')

set.seed(343)
for(i in 1:nsim) {
    require(maxlike)
    cat("sim", i, "\n")
    sim.i <- sim()
    fm.i <- maxlike(~x1 + x2, sim.i$r, sim.i$xy)
    simout[i,] <- coef(fm.i)
}

simout

op <- par(mfrow=c(3,1))
hist(simout[,1])
abline(v=-1, lwd=2, col=4)
hist(simout[,2])
abline(v=1, lwd=2, col=4)
hist(simout[,3])
abline(v=0.5, lwd=2, col=4)
par(op)








# Demonstration that parameters are not identifiable with only
# a single categorical covariate

sim2 <- function() {
    x1 <- raster(matrix(c(0,1), 100, 100), xmn=0, xmx=100, ymn=0, ymx=100)

    logit.psi <- -1 + 1*x1
    psi <- exp(logit.psi)/(1+exp(logit.psi))

    r <- stack(x1)
    r@layernames <- c("x1")

    pa <- matrix(NA, 100, 100)
    pa[] <- rbinom(1e4, 1, as.matrix(psi))
    pa <- raster(pa, 0, 100, 0, 100)
    xy <- xyFromCell(pa, sample(Which(pa==1, cells=TRUE), 1000))

    return(list(r=r, xy=xy))
}


str(sim2(), 1)


plot(sim2()$r)





nsim2 <- 10
simout2 <- matrix(NA, nsim2, 2)
colnames(simout2) <- c('beta0', 'beta1')

set.seed(343)
for(i in 1:nsim2) {
    require(maxlike)
    cat("sim2", i, "\n")
    sim.i <- sim2()
    fm.i <- maxlike(~x1, sim.i$r, sim.i$xy)
    simout2[i,] <- coef(fm.i)
}

simout2

op <- par(mfrow=c(2,1))
hist(simout2[,1])
abline(v=-1, lwd=2, col=4)
hist(simout2[,2])
abline(v=1, lwd=2, col=4)
par(op)














# cloglog link


sim3 <- function() {
    x1 <- sort(rnorm(100))
    x1 <- raster(outer(x1, x1), xmn=0, xmx=100, ymn=0, ymx=100)

    eta <- -1 + 1*x1
    psi <- 1 - exp(-1 * exp(eta))

    r <- stack(x1)
    r@layernames <- c("x1")

    pa <- matrix(NA, 100, 100)
    pa[] <- rbinom(1e4, 1, as.matrix(psi))
    pa <- raster(pa, 0, 100, 0, 100)
    xy <- xyFromCell(pa, sample(Which(pa==1, cells=TRUE), 1000))

    return(list(r=r, xy=xy))
}


str(sim3(), 1)


plot(sim3()$r)





nsim3 <- 100
simout3 <- matrix(NA, nsim3, 2)
colnames(simout3) <- c('beta0', 'beta1')

set.seed(343)
for(i in 1:nsim3) {
    require(maxlike)
    cat("sim3", i, "\n")
    sim.i <- sim3()
    fm.i <- maxlike(~x1, sim.i$r, sim.i$xy, link="cloglog")
    simout3[i,] <- coef(fm.i)
}

simout3

op <- par(mfrow=c(2,1))
hist(simout3[,1])
abline(v=-1, lwd=2, col=4)
hist(simout3[,2])
abline(v=1, lwd=2, col=4)
par(op)

