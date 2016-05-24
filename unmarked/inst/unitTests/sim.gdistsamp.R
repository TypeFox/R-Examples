

library(unmarked)


sim1 <- function(lambda=5, phi=0.5, shape=20, scale=10, R=100, T=3,
    breaks=seq(0, 50, by=10), survey="pt", detfun="hn")
{
    nb <- length(breaks)
    J <- nb-1
    maxDist <- max(breaks)
    tlength <- 1000
    switch(survey,
           pt = A <- (2*maxDist)^2 / 10000,   # Area (ha) of square
           line = A <- maxDist*2*100 / 10000 # Area (ha) 100m transect
           )
    a <- pi*breaks[2]^2
    for(j in 2:J) {
        a[j] <- pi*breaks[j+1]^2 - sum(a[1:(j-1)])
        }
    u <- a / sum(a)
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda * A) # Individuals within the rectangle
        if(identical(survey, "pt")) {
            X <- runif(M, -maxDist, maxDist)
            Y <- runif(M, -maxDist, maxDist)
            d <- sqrt(X^2+Y^2)
            d <- d[d<=maxDist]
            M <- length(d)
            }
        else if(identical(survey, "line"))
            d <- runif(M, -maxDist, maxDist)
        else
            stop("survey must be either 'point' or 'line'")
        if(length(d) > 0) {
            switch(detfun,
                   hn   = p <- exp(-d^2 / (2 * shape^2)),
                   exp  = p <- exp(-d/shape),
                   haz  = p <- 1-exp(-(d/shape)^-scale),
                   unif = p <- 1
                   )
            cp <- p * phi
        } else next
        for(t in 1:T) {
            # Detection process
            d1 <- d[rbinom(length(d), 1, cp) == 1]
            y[i,,t] <- table(cut(d1, breaks, include.lowest=TRUE))
        }
    }
    y <- matrix(y, nrow=R)
    return(y)
}



sim2 <- function(lambda=5, phi=0.5, shape=20, scale=10, R=100, T=3,
    breaks=seq(0, 50, by=10), survey="pt", detfun="hn")
{
    nb <- length(breaks)
    J <- nb-1
    maxDist <- max(breaks)
    tlength <- 1000
    switch(survey,
           pt = A <- pi*maxDist^2 / 10000,   # Area (ha) of circle
           line = A <- maxDist*2*100 / 10000 # Area (ha) 100m transect
           )
    a <- pi*breaks[2]^2
    for(j in 2:J) {
        a[j] <- pi*breaks[j+1]^2 - sum(a[1:(j-1)])
        }
    u <- a / sum(a)
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda * A) # Super-population
        for(t in 1:T) {
            switch(survey,
                pt = {
                    z <- 2*pi*runif(M)
                    u <- runif(M) + runif(M)
                    r <- ifelse(u>1, 2-u, u)
                    X <- maxDist*r*cos(z)
                    Y <- maxDist*r*sin(z)
                    d <- sqrt(X^2+Y^2)
                    d <- d[d<=maxDist]
                    d <- d
                    },
                line = {
                    d <- runif(M, 0, maxDist)
                    })

            # Detection process
            if(length(d) > 0) {
                switch(detfun,
                    hn   = p <- exp(-d^2 / (2 * shape^2)),
                    exp  = p <- exp(-d/shape),
                    haz  = p <- 1-exp(-(d/shape)^-scale),
                    unif = p <- 1
                    )
                cp <- p * phi
                d1 <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d1, breaks, include.lowest=TRUE))
                }
            }
        }
    y <- matrix(y, nrow=R)
    return(y)
}







library(unmarked)

set.seed(3)

breaks <- seq(0, 50, by=10)
T <- 5
umf <- unmarkedFrameGDS(y = sim2(lambda=30, shape=50, phi=0.7,
                                 T=T, breaks=breaks),
                        survey="point", unitsIn="m",
                        dist.breaks=breaks, numPrimary=T)
summary(umf)

system.time(m <- gdistsamp(~1, ~1, ~1, umf, K=200, output="density",
                           starts=c(3, 0, 3))) # 28s

backTransform(m, type="lambda")
backTransform(m, type="phi")
backTransform(m, type="det")



# Point-transect, half-normal
nsim1 <- 2
simout1 <- matrix(NA, nsim1, 3)
colnames(simout1) <- c('lambda', 'phi', 'sigma')
set.seed(4059)
for(i in 1:nsim1) {
    cat("sim1", i, "\n")
    breaks <- seq(0, 50, by=10)
    T <- 5
    y1 <- sim2(lambda=30, shape=50, phi=0.7, R=100, T=T, breaks=breaks)
    umf1 <- unmarkedFrameGDS(y = y1, survey="point",
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m1 <- gdistsamp(~1, ~1, ~1, umf1, output="density", K=100,
                    control=list(trace=TRUE, REPORT=1),
#                    lower=c(-Inf,-5,-Inf), upper=c(Inf, 5, Inf),
                    starts=c(3,0.5,3), se=FALSE)
    e <- coef(m1)
    simout1[i,] <- c(exp(e[1]), plogis(e[2]), exp(e[3]))
    cat("\tbeta.hat =", simout1[i,], "\n")
    }

par(mfrow=c(3, 1))
hist(simout1[,1], xlab=expression(lambda), main="")
abline(v=30, col=4)
hist(simout1[,2], xlab=expression(phi), main=""); abline(v=0.7, col=4)
hist(simout1[,3], xlab=expression(sigma), main=""); abline(v=50, col=4)



# Point-transect, neg exp
nsim2 <- 10
simout2 <- matrix(NA, nsim2, 3)
colnames(simout2) <- c('lambda', 'phi', 'rate')
for(i in 1:nsim2) {
    cat("sim2", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y2 <- sim2(lambda=30, phi=0.7, shape=50, R=100, T=T, breaks=breaks,
               detfun="exp")
    umf2 <- unmarkedFrameGDS(y = y2, survey="point",
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m2 <- gdistsamp(~1,~1,~1, umf2, keyfun="exp", output="density", K=100,
                    starts=c(3,0,3), se=FALSE)
    e <- coef(m2)
    simout2[i,] <- c(exp(e[1]), plogis(e[2]), exp(e[3]))
    cat("\tbeta.hat =", simout2[i,], "\n")
    }

par(mfrow=c(3, 1))
hist(simout2[,1], xlab=expression(lambda), main="")
abline(v=30, col=4, lwd=2)
hist(simout2[,2], xlab=expression(phi), main="")
abline(v=0.7, col=4, lwd=2)
hist(simout2[,3], xlab=expression(rate), main="")
abline(v=50, col=4, lwd=2)




# Point-transect, hazard
nsim <- 10
simout <- matrix(NA, nsim, 4)
colnames(simout) <- c('lambda', 'phi', 'shape', 'scale')
for(i in 1:nsim) {
    cat("sim", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y <- sim2(lambda=30, phi=0.7, shape=30, scale=5, R=100, T=T,
             breaks=breaks, detfun="haz")
    umf <- unmarkedFrameGDS(y = y, survey="point",
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m <- gdistsamp(~1, ~1, ~1, umf, keyfun="hazard", output="density",
                   K=100, starts=c(3, 0, 2, 1), se=FALSE)
    e <- coef(m)
    simout[i,] <- c(exp(e[1]), plogis(e[2]), exp(e[3:4]))
    cat(" beta.hat =", simout[i,], "\n")
    }

par(mfrow=c(3, 1))
hist(simout[,1], xlab=expression(lambda), main="")
abline(v=30, col=4)
hist(simout[,2], xlab=expression(phi), main=""); abline(v=0.7, col=4)
hist(simout[,3], xlab=expression(sigma), main=""); abline(v=30, col=4)




# Point-transect, uniform
nsim4 <- 10
simout4 <- matrix(NA, nsim4, 2)
colnames(simout4) <- c('lambda', 'phi')
for(i in 1:nsim4) {
    cat("sim4", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y4 <- sim2(lambda=20, phi=0.6, R=100, T=T, breaks=breaks,
               detfun="unif",
             survey="pt")
    umf4 <- unmarkedFrameGDS(y = y4, survey="point",
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m4 <- gdistsamp(~1, ~1, ~1, umf4, keyfun="uniform", output="density",
                   K=100, unitsOut="ha", se=FALSE)
    e <- coef(m4)
    simout4[i,] <- c(exp(e[1]), plogis(e[2]))
    cat("\tbeta.hat =", simout4[i,], "\n")
    }

par(mfrow=c(2, 1))
hist(simout4[,1], xlab=expression(lambda), main="")
abline(v=20, col=4)
hist(simout4[,2], xlab=expression(phi), main=""); abline(v=0.6, col=4)










# Line-transect, half-normal
nsim5 <- 10
simout5 <- matrix(NA, nsim5, 3)
colnames(simout5) <- c('lambda', 'phi', 'sigma')
for(i in 1:nsim5) {
    cat("sim5", i, "\n")
    breaks <- seq(0, 50, by=10)
    T <- 5
    y5 <- sim2(lambda=30, phi=0.7, shape=50, R=100, T=T, breaks=breaks,
               survey="line")
    umf5 <- unmarkedFrameGDS(y = y5, survey="line", tlength=rep(100,100),
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m5 <- gdistsamp(~1, ~1, ~1, umf5, output="density", K=100,
                    se=FALSE, starts=c(3, 0, 3))
    e <- coef(m5)
    simout5[i,] <- c(exp(e[1]), plogis(e[2]), exp(e[3]))
    cat("  beta.hat =", simout5[i,], "\n")
    }

par(mfrow=c(3, 1))
hist(simout5[,1], xlab=expression(lambda), main="")
abline(v=30, col=4)
hist(simout5[,2], xlab=expression(phi), main=""); abline(v=0.7, col=4)
hist(simout5[,3], xlab=expression(sigma), main=""); abline(v=50, col=4)



# Line-transect, neg exp
nsim6 <- 10
simout6 <- matrix(NA, nsim6, 3)
colnames(simout6) <- c('lambda', 'phi', 'rate')
for(i in 1:nsim6) {
    cat("sim6", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y6 <- sim2(lambda=30, phi=0.7, shape=30, R=100, T=T, breaks=breaks,
               detfun="exp", survey="line")
    umf6 <- unmarkedFrameGDS(y = y6, survey="line", tlength=rep(100,100),
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m6 <- gdistsamp(~1, ~1, ~1, umf6, keyfun="exp", output="density",
                    K=100, se=FALSE, starts=c(3,0,3))
    e <- coef(m6)
    simout6[i,] <- c(exp(e[1]), plogis(e[2]), exp(e[3]))
    cat("\tbeta.hat =", simout6[i,], "\n")
    }

par(mfrow=c(3, 1))
hist(simout6[,1], xlab=expression(lambda), main="")
abline(v=30, col=4)
hist(simout6[,2], xlab=expression(phi), main=""); abline(v=0.7, col=4)
hist(simout6[,3], xlab=expression(sigma), main=""); abline(v=30, col=4)




# Line-transect, hazard
nsim7 <- 10
simout7 <- matrix(NA, nsim7, 4)
colnames(simout7) <- c('lambda', 'phi', 'shape', 'scale')
for(i in 1:nsim7) {
    cat("sim7", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y7 <- sim2(lambda=30, phi=0.7, shape=30, scale=5, R=100, T=T,
             breaks=breaks, detfun="haz", survey="line")
    umf7 <- unmarkedFrameGDS(y = y7, survey="line", tlength=rep(100,100),
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m7 <- gdistsamp(~1, ~1, ~1, umf7, keyfun="hazard", K=100,
        starts=c(1.5, 0.5, 3, 0), output="density")
    e <- coef(m7)
    simout7[i,] <- c(exp(e[1]), plogis(e[2]), exp(e[3:4]))
    cat("\tbeta.hat =", simout7[i,], "\n")
    }

par(mfrow=c(3, 1))
hist(simout7[,1], xlab=expression(lambda), main="")
abline(v=20, col=4)
hist(simout7[,2], xlab=expression(phi), main=""); abline(v=0.7, col=4)
hist(simout7[,3], xlab=expression(sigma), main=""); abline(v=20, col=4)




# Line-transect, uniform
nsim8 <- 10
simout8 <- matrix(NA, nsim8, 2)
colnames(simout8) <- c('lambda', 'phi')
for(i in 1:nsim8) {
    cat("sim8", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y8 <- sim2(lambda=20, phi=0.7, R=100, T=T, breaks=breaks,
               detfun="unif",
             survey="line")
    umf8 <- unmarkedFrameGDS(y = y8, survey="line", tlength=rep(100,100),
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m8 <- gdistsamp(~1, ~1, ~1, umf8, keyfun="uniform", output="density",
                   K=100, se=FALSE)
    e <- coef(m8)
    simout8[i,] <- c(exp(e[1]), plogis(e[2]))
    }

par(mfrow=c(2, 1))
hist(simout8[,1], xlab=expression(lambda), main="")
abline(v=20, col=4)
hist(simout8[,2], xlab=expression(phi), main=""); abline(v=0.7, col=4)












rcirc <- function(n, R) {
    a <- a2 <- runif(n)
    b <- b2 <- runif(n)
    bla <- b < a
    b2[bla] <- a[bla]
    a2[bla] <- b[bla]
    stopifnot(all(b2 >= a2))
    cbind(b2*R*cos(2*pi*a2/b2), b2*R*sin(2*pi*a2/b2))
}

plot(0,type="n", xlim=c(-55,55), asp=1)
points(rc1 <- rcirc(1000, 50))

rcirc <- function(n, R) {
    t <- 2*pi*runif(n, 0, 1)
    u <- runif(n) + runif(n)
    r <- ifelse(u>1, 2-u, u)
    cbind(R*r*cos(t), R*r*sin(t))
}

plot(0,type="n", xlim=c(-55,55), asp=1)
points(rc1 <- rcirc(1000, 50))



f <- function(r, R=50) 2*r/R^2
u <- runif(10000, 0, 50)
d <- sample(u, 100, replace=TRUE, prob=f(u))
hist(d)
