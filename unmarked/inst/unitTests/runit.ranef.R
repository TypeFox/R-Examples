


# ----------------------------- pcount ----------------------------------

test.ranef.pcount <- function() {

    library(unmarked)
    set.seed(4564)
    R <- 10
    J <- 5
    N <- rpois(R, 3)
    y <- matrix(NA, R, J)
    y[] <- rbinom(R*J, N, 0.5)
    y[1,] <- NA
    y[2,1] <- NA
    K <- 15

    umf <- unmarkedFramePCount(y=y)
    fm <- pcount(~1 ~1, umf, K=K)

    re <- ranef(fm)
    modes <- bup(re, stat="mode")
    CI <- confint(re, level=0.9)
    checkEqualsNumeric(length(modes), R-1)
    checkEqualsNumeric(nrow(CI), R-1)
    checkEqualsNumeric(sum(modes), 42)
    checkEqualsNumeric(colSums(CI), c(33,60))

    df <- as(re, "data.frame")
    ar <- as(re, "array")
    checkEqualsNumeric(nrow(df), 144)
    checkEqualsNumeric(colSums(ar), c(
 0.000000e+00, 0.000000e+00, 8.315871e-01, 1.337063e+00, 1.449232e+00,
 1.804784e+00, 1.839209e+00, 1.137807e+00, 4.478077e-01, 1.226691e-01,
 2.515083e-02, 4.076864e-03, 5.446213e-04, 6.189493e-05, 6.131061e-06,
 5.391007e-07), tolerance=1e-6)

    fm.nb <- update(fm, mix="NB")
    fm.zip <- update(fm, mix="ZIP")

    ar.nb <- as(ranef(fm.nb), "array")
    ar.zip <- as(ranef(fm.zip), "array")

    checkEqualsNumeric(colSums(ar.nb), c(
 0.000000e+00, 0.000000e+00, 8.316904e-01, 1.337039e+00, 1.449088e+00,
 1.804499e+00, 1.839052e+00, 1.137951e+00, 4.480130e-01, 1.227798e-01,
 2.518738e-02, 4.085472e-03, 5.461876e-04, 6.212684e-05, 6.160034e-06,
 5.422348e-07), tolerance=1e-6)

    checkEqualsNumeric(colSums(ar.zip), c(
 0.000000e+00, 0.000000e+00, 8.315404e-01, 1.337022e+00, 1.449214e+00,
 1.804745e+00, 1.839218e+00, 1.137868e+00, 4.478556e-01, 1.226888e-01,
 2.515613e-02, 4.077914e-03, 5.447851e-04, 6.191597e-05, 6.133365e-06,
 5.393215e-07), tolerance=1e-6)


}








# ------------------------------- occu ----------------------------------



test.ranef.occu <- function() {
    set.seed(4564)
    R <- 10
    J <- 5
    z <- rbinom(R, 1, 0.6)
    y <- matrix(NA, R, J)
    y[] <- rbinom(R*J, 1, z*0.7)
    y[1,] <- NA
    y[2,1] <- NA

    x <- y
    x[] <- rnorm(R*J)
    x[3,1] <- NA

    umf <- unmarkedFrameOccu(y=y, obsCovs=list(x=x))
    fm <- occu(~1 ~1, umf)

    re <- ranef(fm)
    modes <- bup(re, stat="mode")
    CI <- confint(re, level=0.95)
    checkEqualsNumeric(length(modes), R-1)
    checkEqualsNumeric(nrow(CI), R-1)
    checkEqualsNumeric(sum(modes), 3)
    checkEqualsNumeric(colSums(CI), c(3,3))

    df <- as(re, "data.frame")
    ar <- as(re, "array")
    checkEqualsNumeric(nrow(df), 18)
    checkEqualsNumeric(colSums(ar), c(5.993957, 3.006043), tolerance=1e-6)

    fmx <- occu(~x ~1, umf)
    arx <- as(ranef(fmx), "array")
    checkEqualsNumeric(colSums(arx), c(5.991553, 3.008447), tolerance=1e-6)

}














# ------------------------------ distsamp -------------------------------



test.distsamp.ranef <- function() {

    set.seed(344)
    lambda <- 10
    sigma <- 20
    npts <- 10
    radius <- 50
    breaks <- seq(0, 50, by=10)
    A <- (2*radius)^2 / 10000 # Area (ha) of square containing circle
    y <- matrix(0, npts, length(breaks)-1)
    N <- integer(npts)
    for(i in 1:npts) {
        M <- rpois(1, lambda * A) # Individuals within the square
        xy <- cbind(x=runif(M, -radius, radius),
                    y=runif(M, -radius, radius))
        d <- apply(xy, 1, function(x) sqrt(x[1]^2 + x[2]^2))
        d <- d[d <= radius]
        N[i] <- length(d)
        if(length(d)) {
            p <- exp(-d^2 / (2 * sigma^2)) # half-normal
            d <- d[rbinom(length(d), 1, p) == 1]
            y[i,] <- table(cut(d, breaks, include.lowest=TRUE))
        }
    }

    umf1 <- unmarkedFrameDS(y = y, survey="point",
                            dist.breaks=breaks, unitsIn="m")
    (m1 <- distsamp(~1 ~1, umf1, starts=c(log(5), log(20))))
    (m2 <- distsamp(~1 ~1, umf1, starts=c(log(5), log(20)),
                    output="abund"))

    re1 <- ranef(m1, K=20)
    re2 <- ranef(m2, K=20)

    checkEquals(mode1 <- bup(re1, stat="mode"), bup(re2, "mode"))
    checkEquals(confint(re1), confint(re2))

    ar1 <- as(re1, "array")

checkEqualsNumeric(colSums(ar1), c(
 0.000000e+00, 2.334960e-01, 8.517322e-01, 1.524261e+00, 1.811577e+00,
 1.691348e+00, 1.421738e+00, 1.085003e+00, 7.119743e-01, 3.898376e-01,
 1.782052e-01, 6.895313e-02, 2.296231e-02, 6.685198e-03, 1.725009e-03,
 3.991224e-04, 8.362689e-05, 1.600128e-05, 2.816112e-06, 4.586885e-07,
 6.951721e-08), tolerance=1e-6)


}











# ------------------------------ multinomPois ----------------------------






test.ranef.multinomPois <- function() {

    # Simulate independent double observer data
    nSites <- 10
    lambda <- 10
    p1 <- 0.5
    p2 <- 0.3
    cp <- c(p1*(1-p2), p2*(1-p1), p1*p2)
    set.seed(9023)
    N <- rpois(nSites, lambda)
    y <- matrix(NA, nSites, 3)
    for(i in 1:nSites) {
        y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
    }

    # Fit model
    observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)
    umf <- unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
                              type="double")
    fm <- multinomPois(~observer-1 ~1, umf)

    # Estimates of random effects
    re <- ranef(fm, K=20)
    ar <- as(re, "array")

checkEqualsNumeric(colSums(ar), c(
 0.0000000000, 0.0837789470, 0.2077360595, 0.3413273644, 0.5043850863,
 0.6810201789, 0.9111517583, 1.2124489471, 1.4341939152, 1.3759825087,
 1.0500832440, 0.7312585749, 0.5381253818, 0.4003005483, 0.2661348133,
 0.1494032715, 0.0705242071, 0.0283773425, 0.0098973528, 0.0030385119,
 0.0008319866), tolerance=1e-6)


}





# ---------------------------- gmultmix ---------------------------------


library(unmarked)
n <- 100  # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods

lam <- 3
phi <- 0.5
p <- 0.3

#set.seed(26)
y <- array(NA, c(n, T, J))
M <- rpois(n, lam)          # Local population size
N <- matrix(NA, n, T)       # Individuals available for detection

for(i in 1:n) {
    N[i,] <- rbinom(T, M[i], phi)
    y[i,,1] <- rbinom(T, N[i,], p)    # Observe some
    Nleft1 <- N[i,] - y[i,,1]         # Remove them
    y[i,,2] <- rbinom(T, Nleft1, p)   # ...
    Nleft2 <- Nleft1 - y[i,,2]
    y[i,,3] <- rbinom(T, Nleft2, p)
    }

y.ijt <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])
umf1 <- unmarkedFrameGMM(y=y.ijt, numPrimary=T, type="removal")

(m1 <- gmultmix(~1, ~1, ~1, data=umf1, K=30))

re <- ranef(m1)
plot(re, layout=c(5,5), xlim=c(-1,20), subset=site%in%1:25)




# ---------------------------- gpcount ---------------------------------


library(unmarked)
test.ranef.gpcount <- function()
{
    y <- matrix(c(1,1,1, 1,0,1, 2,2,2,
                3,2,3, 2,2,2, 1,1,1,
                NA,0,0, 0,0,0, 0,0,0,
                3,3,3, 3,2,3, 2,2,2,
                0,0,0, 0,0,0, 0,0,0), 5, 9, byrow=TRUE)
    siteCovs <- data.frame(x = c(0,2,-1,4,-1))
    obsCovs <- list(o1 = matrix(seq(-3, 3, length=length(y)), 5, 9))
    obsCovs$o1[5,4:6] <- NA
    yrSiteCovs <- list(yr=matrix(c('1','2','2'), 5, 3, byrow=TRUE))
    yrSiteCovs$yr[4,2] <- NA

    umf <- unmarkedFrameGPC(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, numPrimary=3)
    fm <- gpcount(~x, ~1, ~o1, data = umf, K=30)
    re <- ranef(fm)
    checkEqualsNumeric(bup(re, "mode"), c(2,3,0,4,0))

    fm0 <- gpcount(~1, ~1, ~1, data = umf, K=23)
    re0 <- ranef(fm0)
    checkEqualsNumeric(bup(re0, "mode"), c(2,3,0,3,0))
}




# ------------------------------ gdistsamp ------------------------------

test.ranef.gdistsamp <- function() {

    set.seed(36837)
    R <- 10 # number of transects
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
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)
    # Fit the model
    m1 <- gdistsamp(~1, ~1, ~1, umf, output="abund", K=20)
    m2 <- gdistsamp(~1, ~1, ~1, umf, output="density", K=20)

    re1 <- ranef(m1)
    re2 <- ranef(m2)

    ar1 <- as(re1, "array")
    ar2 <- as(re2, "array")

    checkEquals(colSums(ar1), colSums(ar2), tol=1e-5)

    checkEqualsNumeric(colSums(ar1), c(
 0.000000000, 0.000000000, 0.000000000, 0.118002086, 0.307044478,
 0.310241436, 0.239886364, 0.448502098, 0.977448196, 1.436982755,
 1.548430326, 1.359997401, 1.020206300, 0.655381337, 0.372519278,
 0.233418189, 0.218220126, 0.242683228, 0.232040811, 0.174834734,
 0.104160860), tolerance=1e-6)

}

# ----------------------------- colext ----------------------------------




test.ranef.colext <- function() {

    set.seed(7)
    M <- 10
    J <- 3
    T <- 5
    psi <- 0.5
    gamma <- 0.4
    eps <- 0.6
    p <- 0.5
    z <- matrix(NA, M, T)
    y <- array(NA, c(M, J, T))
    z[,1] <- rbinom(M, 1, psi)
    y[,,1] <- rbinom(M*J, 1, z[,1]*p)
    for(t in 1:(T-1)) {
        mu <- ((1-z[,t])*gamma + z[,t]*(1-eps))
        z[,t+1] <- rbinom(M, 1, mu)
        y[,,t+1] <- rbinom(M*J, 1, z[,t+1]*p)
    }

    # Prepare data
    umf <- unmarkedMultFrame(y = matrix(y, M), numPrimary=T)
    summary(umf)

    # Fit model and backtransform
    (m1 <- colext(~1, ~1, ~1, ~1, umf))

    re1 <- ranef(m1)

    plot(re1, xlim=c(-1,2))

    ar1 <- as(re1, "array")


checkEqualsNumeric(colSums(ar1), matrix(c(
 7.587399, 3.711435, 3.138038, 4.677841, 3.121028,
 2.412601, 6.288565, 6.861962, 5.322159, 6.878972), 2, byrow=TRUE),
                   tol=1e-6)



}








# ----------------------------- pcountOpen -------------------------------



test.ranef.pco <- function() {

    set.seed(7)
    M <- 10
    J <- 3
    T <- 5
    lambda <- 5
    gamma <- 0.4
    omega <- 0.6
    p <- 0.5
    N <- matrix(NA, M, T)
    y <- array(NA, c(M, J, T))
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    y[,,1] <- rbinom(M*J, N[,1], p)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma)
        N[,t+1] <- S[,t] + G[,t]
        y[,,t+1] <- rbinom(M*J, N[,t+1], p)
    }

    # Prepare data
    umf <- unmarkedFramePCO(y = matrix(y, M), numPrimary=T)
    summary(umf)

    # Fit model and backtransform
    (m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20))

    re1 <- ranef(m1)
    ar1 <- as(re1, "array")

    #write.table( round(colSums(ar1),6), sep=",", row.names=FALSE,
    #    col.names=FALSE)

    checkEqualsNumeric(colSums(ar1), matrix(c(
0,0,0.819576,1.882332,0.989637,
0.480739,2.960528,3.829933,2.832941,6.712276,
2.288187,2.436992,1.751606,4.066047,2.111822,
2.570878,2.156641,2.933849,1.144065,0.18163,
1.250545,1.555231,0.618936,0.072798,0.00458,
0.971742,0.714645,0.044551,0.001793,5.4e-05,
0.926152,0.156831,0.001519,2.3e-05,0,
0.855516,0.017955,2.9e-05,0,0,
0.469365,0.001134,0,0,0,
0.15106,4.2e-05,0,0,0,
0.030947,1e-06,0,0,0,
0.004375,0,0,0,0,
0.000455,0,0,0,0,
3.7e-05,0,0,0,0,
2e-06,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0), ncol=5, byrow=TRUE), tolerance=1e-6)


    m2 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20, dynamics="trend")
    re2 <- ranef(m2)
    checkEqualsNumeric(bup(re2, "mode")[1,], c(8, 5, 3, 1, 1))

}
