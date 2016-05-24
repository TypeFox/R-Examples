require("tsfa")

if (require("CDNmoney")) {

data("CanadianMoneyData.asof.6Feb2004", package="CDNmoney")

#require("dse")
#require("EvalEst")  # for EstEval

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
             positive.measures=FALSE), Cov=Psi, rng=rngValue10) #, noIC=TRUE)

###  Tests to check that calculated values have not changed       

ML  <- estTSF.ML (simBoblq, 2, BpermuteTarget=Boblq)

if( ! ML$stats$estConverged)       stop("ML estimation did not converge.")
if( ! ML$stats$rotationConverged)  stop("ML rotation did not converge.")

z <- summary(ML)
print(z)

#tst <- c(2.58390678496512294, 4, 0, 1,
#         0.629676597497500046, 8.59534396080311451e-44)
#  early version with fewer stats above

#tst <- c(2.58390678496512294, 4, 0.629676597497500046, 0, 1, 0,
#      0.98904681947754125,  1.02403909700877449,
#      0.958925573040779855,  1.00610635836192164, 1,  1.0033141085616486,
#      0.995955195793569614,  0.978764777916240503,  -5.41609321503487706,
#      110.805499040336599,  93.8054990403365991,  0.170952835443762247,
#      0.176325538500333667)
#  early version with more stats and different order above
tst <- c(2.58390678496512294, 4,  0.629676597497500157, 0, 0, 1, 1,   
       1.0033141085616486,  
       0.995955195793569725, 0.978764777916241058,  -5.41609321503487706,  
       110.805499040336599,  93.8054990403365991,  0.170952835443762247,  
       0.176325538500333667)

 if( fuzz < max(abs(z$fitStats - tst ))) {
    cat("Calculated value is not the same as test value. Value:\n")
    printTestValue(z$fitStats, digits=18)
    cat("difference:\n")
    print(z$fitStats - tst, digits=18)
    all.ok <- FALSE  
    } 


# using early version of GPA (which had maxit=500, eps=1e-5, and possibly did
#   normalizing by default)
# tst <- t(matrix(c( 
#       9.3028072143660356,   0.9176845844525813,
#      17.7502527518722424, -27.3868246616804889,
#       4.2782904551983387,  -3.9956142796066030,
#      41.5374289408217194,  10.086789, #10.0867884770525524
#       6.5686841941148506,  21.8380378686635943,
#      10.3887015530108382,  41.9013875186998987), 2,6))
      
# tst <- t(matrix(c( 
#       9.3028155838216442 , 0.917563628847356649, 
#      17.749882414386537  , -27.387067535568125 , 
#       4.27823573538771651, -3.99567173827507194, 
#      41.5375457163001727 , 10.0862508890412137 , 
#       6.56897103994587184, 21.8379612502913183 , 
#      10.389252837922518  , 41.9012693967163017),  2,6))
#  above was with noIC=TRUE
 tst <- t(matrix(c( 
      8.51204715949594437 , 5.29270407581362878 , 
      19.4220856186493975 , -1.69458912605677714 , 
      6.67771253717154334 , -3.30147684724299673 , 
      31.2032284494761214 , 24.7525488324103122 , 
      1.02360574296459972 , 28.2016066839309296 , 
      -8.01975878958417709 , 54.2866349204937038),  2,6))

 if( fuzz < max(abs(loadings(ML) - tst ))) {
    cat("Calculated value is not the same as test value. Value:\n")
    printTestValue(loadings(ML), digits=18)
    cat("difference:\n")
    print(loadings(ML) - tst, digits=18)
    all.ok <- FALSE  
    } 

# tst <- diag(c(79.01113369909211, 746.40304214354558, 100.13891799114583,
#              395.41379672414382, 3409.11389926934453, 15808.93911608403869))
#  above was with noIC=TRUE
 tst <- diag(c(60.3759759290521671,  1526.00374502951945,  80.3896977488930844,
             541.907521365473826,  3295.48462889734583,  13784.0023152529884))

 if( fuzz < max(abs( TSFmodel(ML)$Omega - tst ))) {
    cat("Calculated value is not the same as test value. Value:\n")
    printTestValue(diag(TSFmodel(ML)$Omega), digits=18)
    cat("difference:\n")
    print(TSFmodel(ML)$Omega - tst, digits=18)
    all.ok <- FALSE  
    }

       
# tst <- t(matrix(c(
#      0.019207692269601, 0.0054898845414494, 0.0087895404088425,
#          0.016512272986346, 5.792619859976e-05, 5.421303133591e-07,
#      0.00216155825722, -0.02315302291361, -0.025873534916491,
#          0.011081480506391, 0.003791393906663, 0.0015745247861089), 6, 2))
#  above was with noIC=TRUE
 tst <- t(matrix(c(
      0.0298662762879307632 , 0.00611800907206768753 , 0.0527539544878877567,
       0.00845396103861114842, -0.00310585488073448438, -0.00175378114817218724, 
      0.0147795889144129913 , -0.00569038907180073866 , -0.0634655440408990468,
       0.0137221905746111837 , 0.00654363825436348953 , 0.00328631180520061548), 6, 2))

 if( fuzz < max(abs(ML$LB - tst ))) {
    cat("Calculated value is not the same as test value. Value:\n")
    printTestValue(ML$LB, digits=18)
    cat("difference:\n")
    print(ML$LB - tst, digits=18)
    all.ok <- FALSE  
    } 


cat("tests completed.\n")


if (! all.ok) stop("some tests FAILED")
} else cat("CDNmoney data not available. Tests skipped.\n")
