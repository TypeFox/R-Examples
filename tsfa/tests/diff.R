# make sure diff(factors(x))  == factors(diff(x))
#  and     diff(explained(x)) == explained(diff(x))

require("tsfa")

if (require("CDNmoney")) {

data("CanadianMoneyData.asof.6Feb2004", package="CDNmoney")

fuzz <- 1e-6 
all.ok <- TRUE  

### Construct data

cpi <- 100 * M1total / M1real
seriesNames(cpi) <- "CPI"
popm <- M1total / M1PerCapita
seriesNames(popm) <- "Population of Canada"

z <- tframed(tbind(
    MB2001,
    MB486 + MB452 + MB453 ,
    NonbankCheq,
    MB472 + MB473 + MB487p,
    MB475,
    NonbankNonCheq + MB454 + NonbankTerm + MB2046 + MB2047 + MB2048 +
    MB2057 + MB2058 + MB482),
    names=c("currency", "personal cheq.", "NonbankCheq",
    "N-P demand & notice", "N-P term", "Investment")
    )

z <- tfwindow(z, start=c(1986,1))
if( all(c(2003,12) ==end(z))) z <-tfwindow(z, end=c(2003,11))
MBcomponents <- 1e8 * z/matrix(tfwindow(popm * cpi,tf=tframe(z)),Tobs(z),6)

### Specify "true" parameters and factors

Omega <- diag(c(72.633490218431234, 1233.026245431895177, 87.337721037020572,
              629.392699084312198, 3967.982989812266169, 12163.258995566555313))

Boblq <- t(matrix(c( 
     8.8424730199723260,   5.2034757439511159,
    23.8239003553122046, -12.5767326858555819,
     5.1878379834837549,  -1.9710231572687940,
    36.7834439249370746,  16.9430526918934632,
    -2.8494845070603847,  31.0224248853059343,
     2.6047417719514878,  47.6304267232332990), 2,6))

PhiOblq <- t(matrix(c(
    1.0000000000000002220, 0.0094910545788177599,
    0.0094910545788177599, 1.0000000000000002220),2,2))

etaBart <- MBcomponents %*% solve(Omega) %*% Boblq %*% (
              solve( t(Boblq) %*% solve(Omega) %*% Boblq ) )

DetaBart <- diff(etaBart, lag=1)
SDE      <- cov(DetaBart)       
RR1 <- chol(SDE)      # upper triangular: SDE = RR1' RR1
RR2 <- chol(PhiOblq)  # ditto
PP  <- t(RR2) %*% solve(t(RR1))

etaTrue <- tframed(etaBart %*% t(PP), tf=tframe(MBcomponents))
Psi       <- 0.5 * Omega

rngValue10 <- list(seed=10, kind="Mersenne-Twister", normal.kind="Inversion")

simBoblq  <- simulate(TSFmodel(Boblq, f=etaTrue, 
             positive.measures=FALSE), Cov=Psi, rng=rngValue10, noIC=TRUE)

ML  <- estTSF.ML (simBoblq, 2, BpermuteTarget=Boblq)

 if( fuzz < max(abs(diff(factors(ML)) - factors(diff(ML)) ))) {
    cat("diff(factors(x))  == factors(diff(x))  test failed\n")
    all.ok <- FALSE  
    } 

 if( fuzz < max(abs(tframe(diff(factors(ML))) - tframe(factors(diff(ML))) ))) {
    cat("tframe(diff(factors(x)))  == tframe(factors(diff(x)))  test failed\n")
    all.ok <- FALSE  
    } 


 if( fuzz < max(abs(diff(explained(ML)) - explained(diff(ML)) ))) {
    cat("diff(explained(x))  == explained(diff(x))  test failed\n")
    all.ok <- FALSE  
    } 

 if( fuzz < max(abs(tframe(diff(explained(ML))) - tframe(explained(diff(ML))) ))) {
    cat("tframe(diff(explained(x)))  == tframe(explained(diff(x)))  test failed\n")
    all.ok <- FALSE  
    } 


cat("tests completed.\n")


if (! all.ok) stop("some tests FAILED")
} else cat("CDNmoney data not available. Tests skipped.\n")
