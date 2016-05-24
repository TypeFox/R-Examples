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
             positive.measures=FALSE), Cov=Psi, rng=rngValue10)#, noIC=TRUE)

###  Tests to check that calculated values have not changed

#QML <- estTSF.QML(simBoblq, 2, BpermuteTarget=Boblq)
QML <- estTSF.ML(simBoblq, 2, BpermuteTarget=Boblq, normalize=FALSE)
       
 tst <- t(matrix(c( 
      8.84048060500597721, -1.22368365056050843,
      1.23157276856511477, -32.6686251750833492,
      1.69497918731446817, -5.16173483155516344,
       42.687216287273408, 0.771317526734934256,
       17.613686024267615, 21.1957659664023339 ,
      31.8084447635782901, 41.1878334274834899), 2,6))

# on breman I get (which is 2.8e-4 from above)
 tst <- t(matrix(c(
      8.84048474869231171 , -1.22369174745068032 , 
      1.23177754092009217 , -32.6685583644094066 , 
      1.69501094026290344 , -5.16172608363218544 , 
      42.6871943231400905 , 0.771264465211652483 , 
      17.6135457821184644 , 21.1957004268134739 , 
      31.8081732134605595 ,  41.187708985369504), 2,6))
#  above was with noIC=TRUE
 tst <- t(matrix(c(
      10.1981977120380662 , -0.975592050780101139 , 
      11.9623533659552432 , -12.8526682858124488 , 
      1.88053144050057131 , -6.51704920005049004 , 
      41.7775952479288364 , 0.552653271591319717 , 
      23.8594010938091756 , 21.1537553833569625 , 
      39.0593458850729291 , 46.6582984817521478), 2,6))

# if( fuzz < max(abs(QML$model$B - tst ))) {
 if( fuzz < max(abs(loadings(QML) - tst ))) {
    cat("Calculated value is not the same as test value in test 8. Value:\n")
    printTestValue(loadings(QML), digits=18)
    cat("difference:\n")
    print(loadings(QML) - tst, digits=18)
    all.ok <- FALSE  
    } 

 tst <- diag(c(79.0111336990921131, 746.403042143545576, 100.138917991145846,
              395.413796724143879,  3409.11389926934453,  15808.9391160840405))
#  above was with noIC=TRUE
 tst <- diag(c(60.3759759290521671,  1526.00374502951945,  80.3896977488930844,
        541.907521365473826,  3295.48462889734583,  13784.0023152529884))

 if( fuzz < max(abs( TSFmodel(QML)$Omega - tst ))) {
    cat("Calculated value is not the same as test value in test 9. Value:\n")
    printTestValue(diag(TSFmodel(QML)$Omega), digits=18)
    cat("difference:\n")
    print(TSFmodel(QML)$Omega - tst, digits=18)
    all.ok <- FALSE  
    }

       
 tst <- t(matrix(c(
  0.0193437516976165125, 0.000275488067858643782, 0.00291522923651014315,
     0.0186675513620792564, 0.000894782856918553476, 0.000348497860818553179,
 -0.00789503495675600644 , -0.0223821887083106262 , -0.0263561813945096812,
     0.00102177524579772109, 0.00318066723415876777, 0.00133280208562447052), 6, 2))
#  above was with noIC=TRUE
 tst <- t(matrix(c(
      0.0312499568411227432, 0.00131405534239187674, 0.00292823273875999351, 
      0.0144119658775881381, 0.00146572563954190552, 0.000589162976325337913,    
      -0.0140980078283146574, -0.0087723553274788469, -0.08536058524511414, 
      0.00244781454332796463, 0.00692031286447671815, 0.00363183819429673026), 6, 2))

 if( fuzz < max(abs(QML$LB - tst ))) {
    cat("Calculated value is not the same as test value in test 8. Value:\n")
    printTestValue(QML$LB, digits=18)
    cat("difference:\n")
    print(QML$LB - tst, digits=18)
    all.ok <- FALSE  
    } 

cat("tests completed.\n")


if (! all.ok) stop("some tests FAILED")
} else cat("CDNmoney data not available. Tests skipped.\n")
