



test.nonparboot.occu <- function() {

    set.seed(3343)
    R <- 20
    J <- 5
    z <- rbinom(R, 1, 0.6)
    y <- matrix(NA, R, J)
    y[] <- rbinom(R*J, 1, z*0.7)
    x1 <- rnorm(R)
    x2 <- y
    x2[] <- rnorm(R*J)
    x2[1,] <- NA
    x2[3,1] <- NA
    umf <- unmarkedFrameOccu(y=y, siteCovs=data.frame(x1=x1),
                             obsCovs=list(x2=x2))
    fm1 <- occu(~1 ~1, umf)
    fm2 <- occu(~x2 ~x1, umf)
    fm1 <- nonparboot(fm1, B=2)
    fm2 <- nonparboot(fm2, B=2)

    checkEqualsNumeric(vcov(fm2, method="nonparboot"), matrix(c(
          0.11708106,  0.07895147, -0.2668053, -0.3152901,
          0.07895147,  0.05323947, -0.1799153, -0.2126101,
         -0.26680534, -0.17991529,  0.6079984,  0.7184859,
         -0.31529012, -0.21261011,  0.7184859,  0.8490516), 4, byrow=TRUE),
          tolerance=1e-5)

}




test.nonparboot.gmultmix <- function() {

    set.seed(34)
    n <- 10       # number of sites
    T <- 3        # number of primary periods
    J <- 3        # number of secondary periods
    lam <- 3
    phi <- 0.5
    p <- 0.3
    #set.seed(26)
    y <- array(NA, c(n, J, T))
    M <- rpois(n, lam)          # Local population size
    N <- matrix(NA, n, T)       # Individuals available for detection
    for(i in 1:n) {
        N[i,] <- rbinom(T, M[i], phi)
        y[i,1,] <- rbinom(T, N[i,], p)    # Observe some
        Nleft1 <- N[i,] - y[i,1,]         # Remove them
        y[i,2,] <- rbinom(T, Nleft1, p)   # ...
        Nleft2 <- Nleft1 - y[i,2,]
        y[i,3,] <- rbinom(T, Nleft2, p)
    }
    y.ijt <- matrix(y, n)
    sc1 <- rnorm(n)
    ysc1 <- matrix(rnorm(n*T), n, T)
    oc1 <- y.ijt
    oc1[] <- rnorm(n*T*J)
    umf1 <- unmarkedFrameGMM(y=y.ijt, siteCovs=data.frame(sc=sc1),
                             yearlySiteCovs=list(ysc=ysc1),
                             obsCovs=list(oc=oc1),
                             numPrimary=T, type="removal")
    fm1 <- gmultmix(~1, ~1, ~1, umf1, K=10)
    fm1 <- nonparboot(fm1, B=2)

    umf2 <- unmarkedFrameGMM(y=y.ijt, siteCovs=data.frame(sc=sc1),
                             obsCovs=list(oc=oc1),
                             numPrimary=T, type="removal")
    fm2 <- gmultmix(~1, ~1, ~1, umf2, K=10)
    fm2 <- nonparboot(fm2, B=2)

    umf3 <- unmarkedFrameGMM(y=y.ijt, siteCovs=data.frame(sc=sc1),
                             numPrimary=T, type="removal")
    fm3 <- gmultmix(~1, ~1, ~1, umf3, K=10)
    fm3 <- nonparboot(fm3, B=2)

    umf4 <- unmarkedFrameGMM(y=y.ijt,
                             numPrimary=T, type="removal")
    fm4 <- gmultmix(~1, ~1, ~1, umf4, K=10)
    fm4 <- nonparboot(fm4, B=2)


}




test.nonparboot.colext <- function() {

    set.seed(343)
    nSites <- 10
    nReps <- 3
    nYears <- 2

    epsilon <- 0.7
    psi <- 0.5
    gamma <- 0.4
    p <- 0.6

    y <- array(NA, c(nSites, nReps, nYears))
    Z <- matrix(NA, nSites, nYears)

    phi <- 1-epsilon

    Z[,1] <- rbinom(nSites, 1, psi)
    for(t in 2:nYears) {
        muZ <- Z[,t-1] * phi + (1 - Z[,t-1]) * gamma
        Z[,t] <- rbinom(nSites, 1, muZ)
    }
    for(j in 1:nReps)
        for(t in 1:nYears)
            y[,j,t] <- rbinom(nSites, 1, Z[,t]*p)

    y <- matrix(y, nSites, nReps*nYears)

    sc <- rnorm(nSites)
    ysc <- matrix(rnorm(nSites*nYears), nSites)
    oc <- matrix(rnorm(nSites*nReps*nYears), nSites)

    umf1 <- unmarkedMultFrame(y=y, siteCovs=data.frame(sc=sc),
                              yearlySiteCovs=list(ysc=ysc),
                              obsCovs=list(oc=oc), numPrimary=nYears)

    m1 <- colext(~1, ~1, ~1, ~1, umf1)
    m1 <- nonparboot(m1, B=2)

    checkEqualsNumeric(vcov(m1, method="nonparboot"),
                       matrix(c(
  1.5313035, -0.4177529,  0.8101477,  0.4844295,
 -0.4177529,  0.1139666, -0.2210153, -0.1321566,
  0.8101477, -0.2210153,  0.4286148,  0.2562911,
  0.4844295, -0.1321566,  0.2562911,  0.1532498), 4, byrow=TRUE),
                       tolerance=1e-6)

    m1.1 <- colext(~sc, ~ysc, ~1, ~oc, umf1)
    m1.1 <- nonparboot(m1.1, B=2)

    sc[1] <- NA
#    ysc[2,1] <- NA
    oc[3,4] <- NA

    umf2 <- unmarkedMultFrame(y=y, siteCovs=data.frame(sc=sc),
                              yearlySiteCovs=list(ysc=ysc),
                              obsCovs=list(oc=oc), numPrimary=nYears)

    m2 <- colext(~sc, ~ysc, ~1, ~oc, umf2)
    m2 <- nonparboot(m2, B=2)



}
