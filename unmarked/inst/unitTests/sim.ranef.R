

# Evaluate bias of BUPs

# ----------------------------- pcount ----------------------------------


library(unmarked)

sim.nmix <- function(R=100, J=5, lambda=5, p=0.7) {
    N <- rpois(R, lambda)
    y <- matrix(NA, R, J)
    y[] <- rbinom(R*J, N, p)
    return(list(y=y, N=N))
}


nsim <- 100
out.nmix <- matrix(NA, nsim, 3)
set.seed(83145)
for(i in 1:nsim) {
    lambda <- 5
    sim.i <- sim.nmix(J=5, lambda=lambda, p=0.5)
    umf <- unmarkedFramePCount(y=sim.i$y)
    K <- 50
    fm <- pcount(~1 ~1, umf, K=K, se=FALSE)
    lam.hat <- exp(coef(fm, type="state"))
    re <- ranef(fm)
    N <- sim.i$N
    N.hat1 <- bup(re, stat="mean")
    N.hat2 <- bup(re, stat="mode")
    bias1 <- mean(N.hat1 - N)
    bias2 <- mean(N.hat2 - N)
    ci <- confint(re)
    cover <- mean(N >= ci[,1] & N <= ci[,2])
    out.nmix[i,] <- c(bias1, bias2, cover)
    cat("sim", i, "\n")
#    cat("  bias =", mean(modes)-lambda, "\n")
}

hist(out.nmix[,1], breaks=20)

colMeans(out.nmix)

plot(re, layout=c(5,5))

plot(N.hat1, N); abline(0,1)
plot(N.hat2, N); abline(0,1)





# ------------------------------- occu ----------------------------------

library(unmarked)




sim.occu <- function(R=100, J=5, psi=0.5, p=0.4) {
    z <- rbinom(R, 1, psi)
    y <- matrix(NA, R, J)
    y[] <- rbinom(R*J, 1, z*p)
    return(list(y=y, z=z))
}


nsim <- 200
out.occu <- matrix(NA, nsim, 3)
set.seed(38845)
for(i in 1:nsim) {
    cat("sim", i, "\n")
    sim.i <- sim.occu(psi=0.8, p=0.4)
    umf <- unmarkedFrameOccu(y=sim.i$y)
    fm <- occu(~1 ~1, umf, se=FALSE)
    re <- ranef(fm)
    z <- sim.i$z
    z.hat1 <- bup(re, stat="mean")
    z.hat2 <- bup(re, stat="mode")
    bias1 <- mean(z.hat1 - z)
    bias2 <- mean(z.hat2 - z)
    ci <- confint(re)
    cover <- mean(z >= ci[,1] & z <= ci[,2])
    out.occu[i,] <- c(bias1, bias2, cover)
}

hist(out.occu[,1])

colMeans(out.occu)

















# ------------------------------ distsamp -------------------------------







sim.ds <- function(lambda=5, shape=20, scale=10, R=100,
    breaks=seq(0, 50, by=10), survey="point", detfun="hn",
    output="density")
{
    nb <- length(breaks)
    J <- nb-1
    maxDist <- max(breaks)
    tlength <- 1000
    if(output=="density") {
        switch(survey,
               point = A <- pi*maxDist^2 / 10000,   # Area (ha) of circle
               line = A <- maxDist*2*100 / 10000 # Area (ha) 100m transect
               )
    } else A <- 1
    a <- pi*breaks[2]^2
    for(j in 2:J) {
        a[j] <- pi*breaks[j+1]^2 - sum(a[1:(j-1)])
    }
    u <- a / sum(a)
    y <- matrix(0, R, J)
    N <- rpois(R, lambda*A)
    for(i in 1:R) {
        switch(survey,
               point = {
                   z <- 2*pi*runif(N[i])
                   u <- runif(N[i]) + runif(N[i])
                   r <- ifelse(u>1, 2-u, u)
                   X <- maxDist*r*cos(z)
                   Y <- maxDist*r*sin(z)
                   d <- sqrt(X^2+Y^2)
                   d <- d[d<=maxDist]
                   d <- d
               },
               line = {
                   d <- runif(N[i], 0, maxDist)
               })

        # Detection process
        if(length(d) > 0) {
            switch(detfun,
                   hn   = p <- exp(-d^2 / (2 * shape^2)),
                   exp  = p <- exp(-d/shape),
                   haz  = p <- 1-exp(-(d/shape)^-scale),
                   unif = p <- 1
                   )
            cp <- p #* phi
            d1 <- d[rbinom(length(d), 1, cp) == 1]
            y[i,] <- table(cut(d1, breaks, include.lowest=TRUE))
        }
    }
    return(list(y=y, N=N))
}





nsim <- 100
out.ds <- matrix(NA, nsim, 3)
set.seed(38845)
for(i in 1:nsim) {
    cat("sim", i, "\n")
    br <- seq(0, 50, by=10)
    sur <- "point"
    ot <- "density"
    lambda <- 20
    sim.i <- sim.ds(lambda=lambda, R=50, breaks=br, survey=sur, output=ot)
    umf <- unmarkedFrameDS(y=sim.i$y, dist.breaks=br, survey=sur,
                           unitsIn="m")
    fm <- distsamp(~1 ~1, umf, output="density", se=FALSE)
    re <- ranef(fm, K=50)
    N <- sim.i$N
    N.hat1 <- bup(re, stat="mean")
    N.hat2 <- bup(re, stat="mode")
    bias1 <- mean(N.hat1 - N)
    bias2 <- mean(N.hat2 - N)
    ci <- confint(re)
    cover <- mean(N >= ci[,1] & N <= ci[,2])
    out.ds[i,] <- c(bias1, bias2, cover)
}

colMeans(out.ds)


plot(N.hat1, N); abline(0,1)
plot(N.hat2, N); abline(0,1)



fmia <- update(fm, output="abund")

re1 <- ranef(fm, K=50)
re2 <- ranef(fmia, K=50)

all.equal(bup(re1), bup(re2), tol=1e-4)



# ------------------------------ multinomPois ----------------------------






# Simulate independent double observer
sim.mn <- function(nSites=50, lambda=10, p1=0.5, p2=0.3) {
    cp <- c(p1*(1-p2), p2*(1-p1), p1*p2)
    N <- rpois(nSites, lambda)
    y <- matrix(NA, nSites, 3)
    for(i in 1:nSites) {
        y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
    }
    return(list(y=y, N=N))
}



nsim <- 500
out.mn <- matrix(NA, nsim, 3)
set.seed(83145)
for(i in 1:nsim) {
    lambda <- 5
    sim.i <- sim.mn(lambda=lambda)
    umf <- unmarkedFrameMPois(y=sim.i$y, type="double")
    fm <- multinomPois(~1 ~1, umf, se=FALSE)
    lam.hat <- exp(coef(fm, type="state"))
    re <- ranef(fm, K=50)
    N <- sim.i$N
    N.hat1 <- bup(re, stat="mean")
    N.hat2 <- bup(re, stat="mode")
    bias1 <- mean(N.hat1 - N)
    bias2 <- mean(N.hat2 - N)
    ci <- confint(re)
    cover <- mean(N >= ci[,1] & N <= ci[,2])
    out.mn[i,] <- c(bias1, bias2, cover)
    if(i %% 10 == 0) cat("sim", i, "\n")
}

hist(out.nmix[,1], breaks=20)

colMeans(out.nmix)

plot(re, layout=c(5,5))

plot(N.hat1, N); abline(0,1)
plot(N.hat2, N); abline(0,1)









# ---------------------------- gmultmix ---------------------------------



sim.gmn <- function(R=50, T, lam=5, phi=0.5, p=0.3) {
    y <- array(NA, c(R, 3, T))
    M <- rpois(R, lam)          # Local population size
    N <- matrix(NA, R, T)       # Individuals available for detection
    for(i in 1:R) {
        N[i,] <- rbinom(T, M[i], phi)
        y[i,1,] <- rbinom(T, N[i,], p)    # Observe some
        Nleft1 <- N[i,] - y[i,1,]         # Remove them
        y[i,2,] <- rbinom(T, Nleft1, p)   # ...
        Nleft2 <- Nleft1 - y[i,2,]
        y[i,3,] <- rbinom(T, Nleft2, p)
    }
    return(list(y=matrix(y,R), M=M))
}






nsim <- 100
out.gmn <- matrix(NA, nsim, 3)
set.seed(831455)
for(i in 1:nsim) {
    R <- 50
    lambda <- 5
    T <- 5
    sim.i <- sim.gmn(R=R, lam=lambda, T=T, p=0.4)
    umf <- unmarkedFrameGMM(y=sim.i$y, numPrimary=T, type="removal")
    fm <- gmultmix(~1, ~1, ~1, umf, se=FALSE, K=40)
    re <- ranef(fm)
    M <- sim.i$M
    M.hat1 <- bup(re, stat="mean")
    M.hat2 <- bup(re, stat="mode")
    bias1 <- mean(M.hat1 - M)
    bias2 <- mean(M.hat2 - M)
    ci <- confint(re)
    cover <- mean(M >= ci[,1] & M <= ci[,2])
    out.gmn[i,] <- c(bias1, bias2, cover)
    if(i %% 1 == 0) {
        cat("sim", i, "\n")
        cat("  lambda =", exp(coef(fm)[1]), "\n")
        cat("  bias1 =", bias1, "\n")
    }
}

hist(out.gmn[,1], breaks=20)

colMeans(out.gmn)

plot(re, layout=c(5,5), subset=site %in% 1:25)

plot(M.hat1, M); abline(0,1)
plot(M.hat2, M); abline(0,1)





umf1 <- unmarkedFrameGMM(y=y.ijt, numPrimary=T, type="removal")

(m1 <- gmultmix(~1, ~1, ~1, data=umf1, K=30))

re <- ranef(m1)
plot(re, layout=c(5,5), xlim=c(-1,20), subset=site%in%1:25)


plot(bup(re, "mode"), M)




# ------------------------------ gdistsamp ------------------------------

set.seed(36837)
R <- 50 # number of transects
T <- 5  # number of replicates
strip.width <- 50
transect.length <- 60 # so that abund != density
breaks <- seq(0, 50, by=10)

lambda <- 10 # Abundance
phi <- 0.6   # Availability
sigma <- 30  # Half-normal shape parameter

J <- length(breaks)-1
y <- array(0, c(R, J, T))
for(i in 1:R) {
    M <- rpois(1, lambda) # Individuals within the 1-ha strip
    for(t in 1:T) {
        # Distances from point
        d <- runif(M, 0, strip.width)
        # Detection process
        if(length(d)) {
            cp <- phi*exp(-d^2 / (2 * sigma^2)) # half-normal w/ g(0)<1
            d <- d[rbinom(length(d), 1, cp) == 1]
            y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
y <- matrix(y, nrow=R) # convert array to matrix

# Organize data
umf <- unmarkedFrameGDS(y = y, survey="line", unitsIn="m",
    dist.breaks=breaks, tlength=rep(transect.length, R), numPrimary=T)
summary(umf)

# Fit the model
m1 <- gdistsamp(~1, ~1, ~1, umf, output="abund", K=50)
summary(m1)
m2 <- gdistsamp(~1, ~1, ~1, umf, output="density", K=50)
summary(m2)


re1 <- ranef(m1)
plot(re1, xlim=c(-1, 30))
re2 <- ranef(m2)
plot(re2, xlim=c(-1, 30))

all.equal(bup(re1), bup(re2), tol=1e-4)
all(confint(re1) == confint(re2))

cbind(bup(re1), bup(re2))


# ----------------------------- colext ----------------------------------



sim.colext <- function(R=50, J=3, T=5, psi=0.5, gamma=0.4, eps=0.6,
                       p=0.5) {
    z <- matrix(NA, R, T)
    y <- array(NA, c(R, J, T))
    z[,1] <- rbinom(R, 1, psi)
    y[,,1] <- rbinom(R*J, 1, z[,1]*p)
    for(t in 1:(T-1)) {
        mu <- ((1-z[,t])*gamma + z[,t]*(1-eps))
        z[,t+1] <- rbinom(R, 1, mu)
        y[,,t+1] <- rbinom(R*J, 1, z[,t+1]*p)
    }
    return(list(y=matrix(y,R), z=z))
}



nsim <- 100
out.colext <- matrix(NA, nsim, 3)
set.seed(83145)
for(i in 1:nsim) {
    R <- 50
    T <- 10
    sim.i <- sim.colext(R=R, T=T)
    umf <- unmarkedMultFrame(y=sim.i$y, numPrimary=T)
    fm <- colext(~1, ~1, ~1, ~1, umf, se=FALSE)
    re <- ranef(fm)
    z <- sim.i$z
    z.hat1 <- bup(re, stat="mean")
    z.hat2 <- bup(re, stat="mode")
    bias1 <- mean(z.hat1 - z)
    bias2 <- mean(z.hat2 - z)
    ci <- confint(re)
    cover <- mean(z >= ci[,1,] & z <= ci[,2,])
    out.colext[i,] <- c(bias1, bias2, cover)
    if(i %% 1 == 0) cat("sim", i, "\n")
}

hist(out.colext[,1], breaks=20)

colMeans(out.colext)

plot(re, layout=c(5,5), xlim=c(-2,3))

plot(z.hat1, z); abline(0,1)
plot(z.hat2, z); abline(0,1)




# ----------------------------- pcountOpen -------------------------------




library(unmarked)
set.seed(7)

sim.pco <- function(R=100, J=3, T=10, lambda=5, gamma=0.4, omega=0.9,
                    p=0.5) {
    N <- matrix(NA, R, T)
    y <- array(NA, c(R, J, T))
    S <- G <- matrix(NA, R, T-1)
    N[,1] <- rpois(R, lambda)
    y[,,1] <- rbinom(R*J, N[,1], p)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(R, N[,t], omega)
        G[,t] <- rpois(R, gamma)
        N[,t+1] <- S[,t] + G[,t]
        y[,,t+1] <- rbinom(R*J, N[,t+1], p)
    }
    return(list(y=matrix(y,R), N=N))
}


nsim <- 10
out.pco <- matrix(NA, nsim, 3)
set.seed(83145)
for(i in 1:nsim) {
    R <- 50
    lambda <- 5
    T <- 10
    sim.i <- sim.pco(R=R, lambda=lambda, T=T)
    umf <- unmarkedFramePCO(y=sim.i$y, numPrimary=T)
    fm <- pcountOpen(~1, ~1, ~1, ~1, umf, se=FALSE, K=20)
    re <- ranef(fm) # really slow
    N <- sim.i$N
    N.hat1 <- bup(re, stat="mean")
    N.hat2 <- bup(re, stat="mode")
    bias1 <- mean(N.hat1 - N)
    bias2 <- mean(N.hat2 - N)
    ci <- confint(re)
    cover <- mean(N >= ci[,1,] & N <= ci[,2,])
    out.pco[i,] <- c(bias1, bias2, cover)
    if(i %% 1 == 0) cat("sim", i, "\n")
}

hist(out.nmix[,1], breaks=20)

colMeans(out.nmix)

plot(re, layout=c(5,5))

plot(N.hat1, N); abline(0,1)
plot(N.hat2, N); abline(0,1)



(fm.nt <- update(fm, dynamics="notrend"))

re <- ranef(fm.nt)

plot(re, layout=c(5,5), subset = site %in% 1:25, xlim=c(-1,10))


