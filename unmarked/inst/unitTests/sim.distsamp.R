



# Simulate with half-normal detection function

simPt <- function(lambda=5, sigma=20, npts=100, radius=50,
    breaks=seq(0, 50, by=10))
{
    A <- (2*radius)^2 / 10000 # Area (ha) of square containing circle
    y <- matrix(0, npts, length(breaks)-1)
    for(i in 1:npts) {
        M <- rpois(1, lambda * A) # Individuals within the square
        # coordinates of each individual
        xy <- cbind(x=runif(M, -radius, radius), y=runif(M, -radius, radius))

        # Distances from each point
        d <- apply(xy, 1, function(x) sqrt(x[1]^2 + x[2]^2))
        d <- d[d <= radius]

        # Detection process
        if(length(d)) {
            p <- exp(-d^2 / (2 * sigma^2)) # half-normal
            d <- d[rbinom(length(d), 1, p) == 1]
            y[i,] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    return(y)
}

colSums(simPt())

set.seed(3)
umf1 <- unmarkedFrameDS(y = simPt(), survey="point",
    dist.breaks=seq(0, 50, by=10), unitsIn="m")
(m1 <- distsamp(~1 ~1, umf1, starts=c(log(5), log(20))))
(m2 <- distsamp(~1 ~1, umf1, starts=c(log(5), log(20)), output="abund"))


checkEqualsNumeric(coef(m1), c(1.813819, 2.893771), tol=1e-5)
checkEquals(exp(coef(m1, type="state")),
    exp(coef(m2, type="state")) / (pi * 50^2 / 10000), tol=0.01)


set.seed(11)
nsims <- 50
simout1 <- matrix(NA, nsims, 2)
lam <- 20
sig <- 30
for(i in 1:nsims) {
    cat("sim", i, "\n")
    umf <- unmarkedFrameDS(y = simPt(lambda=lam, sigma=sig), survey="point",
        dist.breaks=seq(0, 50, by=10), unitsIn="m")
    m <- distsamp(~1 ~1, umf, starts=c(log(lam), log(sig)),
        output="abund")
    simout1[i,] <- exp(coef(m))
    }
hist(simout1[,1]); abline(v=lam*pi*50^2/10000, lwd=2, col=3)
hist(simout1[,2]); abline(v=sig, lwd=2, col=3)






integrate(unmarked:::grhn, 0, 10, sigma=1000)$value * 2 * pi
integrate(unmarked:::grhn, 10, 20, sigma=1000)$value * 2 * pi


fitstats <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2)
    chisq <- sum((observed - expected)^2 / expected)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
    }
pb <- parboot(m1, statistic=fitstats, nsim=200, report=1)










simLine <- function(lambda=5, sigma=20, npts=100,
    breaks=seq(0, 50, by=10))
{
    W <- max(breaks)
    A <- 2*W*100/10000 # Area (ha) of rectangle containing 100m line
    y <- matrix(0, npts, length(breaks)-1)
    for(i in 1:npts) {
        N <- rpois(1, lambda * A) # Individuals within the square
        # distance from the line
        d <- runif(N, 0, W)

        # Detection process
        if(length(d) > 0) {
            p <- exp(-d^2 / (2 * sigma^2)) # half-normal
            d <- d[rbinom(length(d), 1, p) == 1]
            y[i,] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    return(y)
}

simLine()




set.seed(7)
nsims <- 100
simout2 <- matrix(NA, nsims, 2)
lam <- 20
sig <- 30
for(i in 1:nsims) {
    cat("sim", i, "\n"); flush.console()
    y.sim <- simLine(lambda=lam, sigma=sig)
    umf <- unmarkedFrameDS(y = y.sim, survey="line",
        dist.breaks=seq(0, 50, by=10), unitsIn="m",
        tlength=rep(100, nrow(y.sim)))
    m <- distsamp(~1 ~1, umf, starts=c(log(lam), log(sig)), rel.tol=1e-3)
    simout2[i,] <- exp(coef(m))
    }
hist(simout2[,1]); abline(v=lam, lwd=2, col=3)
hist(simout2[,2]); abline(v=sig, lwd=2, col=3)













# Integrate fails if function is flat over most of its range


grhaz <- unmarked:::grhaz
curve(grhaz(x, shape=100, scale=10), 0, 50)
curve(grhaz(x, shape=1, scale=1), 0, 50)

str(integrate(grhaz, 0, 50, shape=1, scale=1, abs.tol=1e-4, stop.on.error=F))
str(integrate(grhaz, 0, Inf, shape=1, scale=1, abs.tol=1e-1, stop.on.error=F))


curve(grhaz(x, shape=0.1, scale=1), 400, 450)

str(integrate(grhaz, 400, 450, shape=0.1, scale=1, abs.tol=1e-4,
    stop.on.error=F))
str(integrate(grhaz, 300, Inf, shape=0.1, scale=1, abs.tol=1e-1,
    stop.on.error=F))



str(integrate(grhaz, 400, 450, shape=10, scale=10, abs.tol=1e-4,
    stop.on.error=F))


curve(gxhaz(x, shape=5, scale=1), 0, 50)
integrate(gxhaz, 0, 20, shape=5, scale=1)$value

increment <- 2
sum(gxhaz(seq(0, 20, by=increment), shape=5, scale=1) * increment)

increment <- 1
sum(gxhaz(seq(0, 20, by=increment), shape=5, scale=1) * increment)

increment <- 0.2
sum(gxhaz(seq(0, 20, by=increment), shape=5, scale=1) * increment)

increment <- 0.001
sum(gxhaz(seq(0, 20, by=increment), shape=5, scale=1) * increment)









# numeric vs analytic integration

# ------------------------------ Point transects ----------------------

# Half-normal

sigma <- 30
br <- seq(0, 50, 10)

# numeric
hn1 <- integrate(grhn, br[1], br[2], sigma=sigma)$value
a1 <- pi*br[2]^2
p1 <- 2*pi*hn1/(a1)
for(i in 2:(length(br)-1)) {
    hn1[i] <- integrate(grhn, br[i], br[i+1], sigma=sigma)$value
    a1[i] <- (pi*br[i+1]^2) - sum(a1[1:(i-1)])
    p1[i] <- 2*pi*hn1[i]/a1[i]
}
p1


# analytical
hn2 <- sigma^2*(1 - exp(-br[2]^2/(2*sigma^2))) -
    sigma^2*(1 - exp(-br[1]^2/(2*sigma^2)))
a2 <- pi*br[2]^2
p2 <- 2*pi*hn2/(a2)
for(i in 2:(length(br)-1)) {
    hn2[i] <- integrate(grhn, br[i], br[i+1], sigma=sigma)$value
    a2[i] <- (pi*br[i+1]^2) - sum(a2[1:(i-1)])
    p2[i] <- 2*pi*hn2[i]/a2[i]
}
p2

p1






# Negative exponential

rate <- 30
br <- seq(0, 50, 10)
a <- pi*br[2]^2

# numeric
nexp1 <- integrate(unmarked:::grexp, br[1], br[2], rate=rate)$value
p1 <- 2*pi*nexp1/a
for(i in 2:(length(br)-1)) {
    nexp1[i] <- integrate(unmarked:::grexp, br[i], br[i+1],
                          rate=rate)$value
    a[i] <- (pi*br[i+1]^2) - sum(a[1:(i-1)])
    p1[i] <- 2*pi*nexp1[i]/a[i]
}
p1


# analytical

f <- function(x, a) exp(x*a)*x
f(10, a=-1/30)
grexp(10, 30)

f(-10, a=1/30)
grexp(10, 30)

fp <- function(x, a) (x/a - 1/a^2)*exp(x*a)
fp(10, a=-1/30)

fp(10, a=.5)
integrate(f, 0, 10, a=0.5)$value

fp(10, a=-.5)
integrate(f, 0, 10, a=-0.5)$value
fp2 <- function(x, a) (a*x - 1)*exp(a*x) / a^2
fp2(10, a=-.5)


# Integral of exp(x*a)*x
# f=exp(x*a)
# F=x*(1-exp(x*a))
# g=x
# gp=x^2/2
# rate*(1 - exp(-br[2]/rate))
fp <- function(x, a) {
    g <- x
    gp <- x^2/2
    F <- x*(1-exp(x*a))
    F*g - 0.5*3*x^2#F*g - \int F*gp
}


fp2 <- function(x, a) (x/abs(a) - 1/a^2)*exp(x*abs(a))
fp2(10, a=-.5)

exp(r*(-1/rate))*r

(10/(1/rate) - 1/(1/rate)^2)*exp(10*(-1/rate))

(nexp2 <- rate^2*(1 - exp(-br[2]^2/(rate))) -
    rate^2*(1 - exp(-br[1]^2/(rate))))
p2 <- 2*pi*nexp2/a
for(i in 2:(length(br)-1)) {
    nexp2[i] <- integrate(grhn, br[i], br[i+1], rate=rate)$value
    p2[i] <- 2*pi*nexp2[i]/a[i]
}
p2

p1

nexp1
nexp2


30^2*(1-exp(-10/30)) - 0.5*3*30^2*(1-exp(-10/30))


# ------------------------------ Line transects ----------------------



# Half-normal

sigma <- 30
br <- seq(0, 50, 10)
a <- br[2]

# numeric
hn1 <- integrate(gxhn, br[1], br[2], sigma=sigma)$value
p1 <- hn1/a
for(i in 2:(length(br)-1)) {
    hn1[i] <- integrate(gxhn, br[i], br[i+1], sigma=sigma)$value
    a[i] <- br[i+1] - sum(a[1:(i-1)])
    p1[i] <- hn1[i]/a[i]
}
p1


# analytical
hn2 <- sigma^2*(1 - exp(-br[2]^2/(2*sigma^2))) -
    sigma^2*(1 - exp(-br[1]^2/(2*sigma^2)))
p2 <- hn2/a
for(i in 2:(length(br)-1)) {
    hn2[i] <- sigma*(1 - exp(-br[i+1]/sigma)) -
    sigma*(1 - exp(-br[i]/sigma))
    p2[i] <- hn2[i]/a
}
p2

p1

hn1
hn2





# Negative exponential

rate <- 30
br <- seq(0, 50, 10)
a <- br[2]

# numeric
nexp1 <- integrate(gxexp, br[1], br[2], rate=rate)$value
p1 <- nexp1/a
for(i in 2:(length(br)-1)) {
    nexp1[i] <- integrate(gxexp, br[i], br[i+1], rate=rate)$value
    a[i] <- br[i+1] - sum(a[1:(i-1)])
    p1[i] <- nexp1[i]/a[i]
}
p1


# analytical
nexp2 <- rate*(1 - exp(-br[2]/rate)) -
    rate*(1 - exp(-br[1]/rate))
p2 <- nexp2/a
for(i in 2:(length(br)-1)) {
    nexp2[i] <- rate*(1 - exp(-br[i+1]/rate)) -
    rate*(1 - exp(-br[i]/rate))
    p2[i] <- nexp2[i]/a
}
p2

p1

nexp1
nexp2
