#  Tests of calculated "true values" for simulations.
#  See notes in ~/papers/InProgress/TimeSeriesFA/calculations.R

require("tsfa")

if (require("CDNmoney")) {

data("CanadianMoneyData.asof.6Feb2004", package="CDNmoney")

# for monte carlo
require("setRNG")
require("dse")
require("EvalEst")  # for EstEval

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




##poppar

DX <- diff(MBcomponents)
stds <- sqrt(diag(cov(DX)))

poppar <- factanal(factors = 2, covmat = cov(DX), n.obs = Tobs(DX),
            scores = "none", rotation = "none")

# Quartimin (= oblimin with parameter 0) rotation with Kaiser normalization
rownorms  <- sqrt(rowSums(poppar$loadings^2))

KBorth    <- diag(1/rownorms) %*% poppar$loadings
rotpoppar <- GPFoblq(A = KBorth, Tmat = diag(1, 2), method="quartimin")
Boblq     <- diag(stds * rownorms) %*% rotpoppar$loadings
PhiOblq   <- rotpoppar$Phi
Omega     <- diag(stds * poppar$uniquenesses * stds)
Psi       <- 0.5 * Omega

# So the true parameters in the simulation are
  Boblq
  PhiOblq
  Omega
  Psi

SigmaOblq <- Boblq %*% PhiOblq %*% t(Boblq) + Omega

etaBart <- MBcomponents %*% solve(Omega) %*% Boblq %*% (
              solve( t(Boblq) %*% solve(Omega) %*% Boblq ) )

DetaBart <- diff(etaBart, lag=1)
SDE      <- cov(DetaBart)       
RR1 <- chol(SDE)      # upper triangular: SDE = RR1' RR1
RR2 <- chol(PhiOblq)  # ditto
PP  <- t(RR2) %*% solve(t(RR1))

etaTrue <- tframed(etaBart %*% t(PP), tf=tframe(MBcomponents))


###  Tests to check that calculated values have not changed

#  NB  fuzz is relaxed even more further below
fuzz <- 1e-5 # this is pretty relaxed compared to my usual tests, eventually
#              tighten up with values printed to higher precision.

# ######################################################################
# # How well do the means fit?
 
 meanDX  <- colMeans(diff(MBcomponents))
 SigmaOblq <- Boblq %*% PhiOblq %*% t(Boblq) + Omega

# kappaOblq <- solve(t(Boblq) %*% solve(SigmaOblq) %*% Boblq) %*%
#		      t(Boblq) %*% solve(SigmaOblq) %*% t(meanDX)
 kappaOblq <- solve(t(Boblq) %*% solve(SigmaOblq, Boblq),
		      t(Boblq)) %*% solve(SigmaOblq) %*% meanDX
 
 meanhat <- t(Boblq %*% kappaOblq)
 

######################################################################

 tst <- c(  2.882879, -0.3320896, 0.1614583, 10.56027, 9.028839, 14.88725)
 if( fuzz < max(abs(meanhat - tst))) {
    cat("Calculated value is not the same as test value in test 1. Value:\n")
    printTestValue(meanhat, digits=18)
    cat("difference:\n")
    print(meanhat - tst, digits=18)
    all.ok <- FALSE  
    }
    
 tst <- c(  1.623936, 0.8437919, 2.553238, 10.21236, 9.243064, 56.72988)
 if( fuzz < max(abs( meanDX - tst ))) {
    cat("Calculated value is not the same as test value in test 2. Value:\n")
    printTestValue(meanDX, digits=18)
    cat("difference:\n")
    print(meanDX - tst, digits=18)
    all.ok <- FALSE  
    }


 tst <- c(  77.524138564771, -139.35681520788, -93.676330196783, 3.4067233632789,
      -2.3176922549429, -73.757662639795)
 if( fuzz < max(abs( 100*(meanhat-meanDX)/meanDX - tst ))) {
    cat("Calculated value is not the same as test value in test 3. Value:\n")
    printTestValue(100*(meanhat-meanDX)/meanDX, digits=18)
    cat("difference:\n")
    print(100*(meanhat-meanDX)/meanDX - tst, digits=18)
    all.ok <- FALSE  
    }

 
              
# Check correctness
DetaTrue <- diff(etaTrue, lag=1)
SigDeta  <- cov(DetaTrue)
CheckFax <- 100 * (SigDeta - PhiOblq)/PhiOblq   # relative difference in %

 if( fuzz < max(abs( CheckFax ))) {
    cat("Calculated value is not the same as test value in test 4. Value:\n")
    printTestValue(CheckFax, digits=18)
    cat("difference:\n")
    print(CheckFax - tst, digits=18)
    all.ok <- FALSE  
    }



 tst <- c(  0.40629016392777, 0.63131919291667, 0.74051211937755, 0.27589064829569,
         0.80375582540918, 0.84226727885335)
 if( fuzz < max(abs(  poppar$uniquenesses - tst))) {
    cat("Calculated value is not the same as test value in test 5. Value:\n")
    printTestValue(poppar$uniquenesses, digits=18)
    cat("difference:\n")
    print(poppar$uniquenesses - tst, digits=18)
    all.ok <- FALSE  
    }


 tst <- t(matrix(c( 
       0.76269855837340,  0.109545890987616,
       0.38940769411444, -0.465876126560372,
       0.37235218322458, -0.347621664499332,
       0.85015261178239,  0.036737207574017,
       0.13286428348581,  0.422603915908731,
       0.17306303686227,  0.357462518428758), 2,6))
 if( fuzz < max(abs(  poppar$loadings - tst ))) {
    cat("Calculated value is not the same as test value in test 6. Value:\n")
    printTestValue(poppar$loadings[,], digits=18)
    cat("difference:\n")
    print(poppar$loadings[,] - tst, digits=18)
    all.ok <- FALSE  
    }


#poppar
#		Factor1 Factor2
#SS loadings	  1.642   0.658
#Proportion Var   0.274   0.110
#Cumulative Var   0.274   0.383

fuzz <- 1e-6 

## simulation values

#tst <- diag(c(72.63349, 1233.026, 87.33772, 629.3927, 3967.983, 12163.26))
tst <- diag(c(72.633490218431234, 1233.026245431895177, 87.337721037020572,
              629.392699084312198, 3967.982989812266169, 12163.258995566555313))
 if( fuzz < max(abs( Omega - tst ))) {
    cat("Calculated value is not the same as test value in test 7. Value:\n")
    printTestValue(diag(Omega), digits=18)
    cat("difference:\n")
    print(Omega - tst, digits=18)
    all.ok <- FALSE  
    }


# using early version of GPA 
# tst <- t(matrix(c( 
#	8.8424730937083,   5.2034756517999,
#      23.8239001815301, -12.5767329305997,
#	5.1878379563914,  -1.9710232106657,
#      36.7834441658135,  16.9430523091825,
#      -2.8494840729111,  31.0224249106556,
#	2.6047424392098,  47.6304266899614), 2,6))

 tst <- t(matrix(c( 
      8.84252375520804534,   5.20339025247098252, 
     23.8237818511178361,  -12.576960598600591, 
      5.18781954415890123,  -1.97107285697162915, 
     36.7836098529645881,   16.9426974878479477, 
     -2.8491870985089478,   31.0224494555067913, 
      2.60519903617723481,  47.6303973396630411), 2,6))

 if( fuzz < max(abs( Boblq - tst ))) {
    cat("Calculated value is not the same as test value in test 7.5. Value:\n")
    printTestValue(Boblq, digits=18)
    cat("difference:\n")
    print(Boblq - tst, digits=18)
    all.ok <- FALSE  
    } 


cat("tests completed.\n")

if (! all.ok) stop("some tests FAILED")
} else cat("CDNmoney data not available. Tests skipped.\n")
