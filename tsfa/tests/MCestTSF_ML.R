
#  See notes in file calculations.R for TimeSeriesFA paper 

require("tsfa")

if (require("CDNmoney")) {

data("CanadianMoneyData.asof.6Feb2004", package="CDNmoney")

# for monte carlo
require("dse")
require("EvalEst")  # for EstEval

fuzz <- 1e-8
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
Psi       <- 0.5 * Omega

etaTrue <- tframed(etaBart %*% t(PP), tf=tframe(MBcomponents))


rngValue10 <- list(seed=10, kind="Mersenne-Twister", normal.kind="Inversion")

cat("Starting the trucated monte carlo experiment test...\n")

EE.ML5 <- EstEval(TSFmodel(Boblq, f=etaTrue, positive.measures=FALSE),
    replications=5, rng=rngValue10, quiet=FALSE,
    simulation.args=list(Cov=Psi),
    estimation="estTSF.ML", 
    estimation.args=list(2, BpermuteTarget=Boblq),
    criterion ="TSFmodel")

mc5SummaryML <- summary(EE.ML5)

# the next were not tested with early version of GPA  

 tst <- c(-6.35811456287209165,  0.0607392955938533957)
#  above was with noIC=TRUE
 tst <- c(-2.92588953587939571,  -4.71404513263362457)

 if( fuzz < max(abs( (mc5SummaryML$meanhatf.error)  - tst ))) {
    cat("Calculated value is not the same as test value in test 19. Value:\n")
    printTestValue(mc5SummaryML$meanhatf.error, digits=18)
    cat("difference:\n")
    print(mc5SummaryML$meanhatf.error - tst, digits=18)
    all.ok <- FALSE  
    }

 tst <-  c(11.1176063210974636,  18.1248814664440161)
#  above was with noIC=TRUE
 tst <-  c(8.60111941828624538,  24.7872684764427618)

 if( fuzz < max(abs( mc5SummaryML$meanSDhatf  - tst ))) {
    cat("Calculated value is not the same as test value in test 21. Value:\n")
    printTestValue(mc5SummaryML$meanSDhatf, digits=18)
    cat("difference:\n")
    print(mc5SummaryML$meanSDhatf - tst, digits=18)
    all.ok <- FALSE  
    }

 tst <- t(matrix(c( 
      -1.07625157111048164 , 0.703342731082228667, 
      -1.2178110662245345  , 4.70567370637503757 , 
      -0.167039105282395184, 0.226947650291669634, 
      -6.57680144713607007 , 6.63837150437990431 , 
      -7.64998746817823871 , -1.11899415928688128 , 
      -4.34709518945724671 , 3.75421190419606887 ), 2,6))
#  above was with noIC=TRUE
 tst <- t(matrix(c( 
      -1.06570250190681648 , -0.0130131676780012029 , 
      -1.85077437881246709 , 3.15697064361194002 , 
      -0.421861915762349327 , -0.337332002481196769 , 
      -2.95530986474609847 , -0.00875023040947908726 , 
      1.05168217191284175 , 2.10624959737685558 , 
      2.07817788622315858 , -7.99355084897341328), 2,6))

 if( fuzz < max(abs( mc5SummaryML$meanhatB.error - tst ))) {
    cat("Calculated value is not the same as test value in test 23. Value:\n")
    printTestValue(mc5SummaryML$meanhatB.error, digits=18)
    cat("difference:\n")
    print(mc5SummaryML$meanhatB.error - tst, digits=18)
    all.ok <- FALSE  
    }

 tst <- t(matrix(c( 
      2.09171817126895432 , 3.10870214359533614 , 
      3.31424639812060606 , 9.86511334745889101 , 
      0.510677029668560656 , 1.93445199388070876 , 
      8.13415536978554243 , 11.7227922878270903 , 
     13.0294442536303503 , 6.96160177068301955 , 
     10.7630330480164904 , 19.6839038208302632), 2,6))
#  above was with noIC=TRUE
 tst <- t(matrix(c( 
      2.00786152578053922 , 3.51641487208566783 , 
      5.55688084114462377 , 9.53701794115509216 , 
       1.2949230299373764 , 1.46529585264552353 , 
      5.85468954626221016 , 10.3962361962235175 , 
        14.68586418448686 , 7.04445496431349572 , 
      11.7090617911295691 , 19.4653799970791646), 2,6))

 if( fuzz < max(abs( mc5SummaryML$SDhatB - tst ))) {
    cat("Calculated value is not the same as test value in test 24. Value:\n")
    printTestValue(mc5SummaryML$SDhatB, digits=18)
    cat("difference:\n")
    print(mc5SummaryML$SDhatB - tst, digits=18)
    all.ok <- FALSE  
    }


 fuzz <- 1e-5
 
 tst <- c(-0.0268800596792208475,  -0.0446045635023343962)
#  above was with noIC=TRUE
 tst <- c( -0.00762399780916277366,  -0.0437505013835693748)

 if( fuzz < max(abs(  (mc5SummaryML$meanhatDf.error)  - tst ))) {
    cat("Calculated value is not the same as test value in test 27. Value:\n")
    printTestValue(mc5SummaryML$meanhatDf.error, digits=18)
    cat("difference:\n")
    print(mc5SummaryML$meanhatDf.error - tst, digits=18)
    all.ok <- FALSE  
    }

 tst <- c( 0.689731516283344948,  0.80395876793842791)
#  above was with noIC=TRUE
 tst <- c(0.622493027392810316,  0.878995874209845596)

 if( fuzz < max(abs( mc5SummaryML$meanSDhatDf  - tst ))) {
    cat("Calculated value is not the same as test value in test 29. Value:\n")
    printTestValue(mc5SummaryML$meanSDhatDf, digits=18)
    cat("difference:\n")
    print(mc5SummaryML$meanSDhatDf - tst, digits=18)
    all.ok <- FALSE  
    }



#### this is fairly long 

cat("Starting the main monte carlo experiment...\n")

EE.ML100 <- EstEval(TSFmodel(Boblq, f=etaTrue, positive.measures=FALSE),
    replications=100, rng=rngValue10, quiet=FALSE,
    simulation.args=list(Cov=Psi),
    estimation="estTSF.ML", 
    estimation.args=list(2, BpermuteTarget=Boblq),
    criterion ="TSFmodel")

mc100SummaryML <- summary(EE.ML100)

fuzz <- 1e-5  # there seems to be a fair amount of "wobble" in these results 


# tst <- c(-1.437001,        -7.074382)
# tst <- c(-1.4366371949955, -7.0749007937085)
# using early version of GPA 
# tst <- c(-1.4366366998990,     -7.0749017751486)
  tst <- c(-1.43661535460228751, -7.07498653794059429)
#  above was with noIC=TRUE
  tst <- c(-4.25581051242354835,  -6.55708042849400741) # Linux Mandrake 32bit
  tst <- c(-4.25580918879642311,  -6.55708238950199629) # Solaris
  tst <- c(-4.2558098,            -6.5570813) 

 if( fuzz < max(abs( (mc100SummaryML$meanhatf.error)  - tst ))) {
    cat("Calculated value is not the same as test value in test 19. Value:\n")
    printTestValue(mc100SummaryML$meanhatf.error, digits=18)
    cat("difference:\n")
    print(mc100SummaryML$meanhatf.error - tst, digits=18)
    all.ok <- FALSE  
    }

# tst <- c(10.72659,        20.42743)
# using early version of GPA 
# tst <-  c(10.726408029255,     20.427522660190)
  tst <-  c(10.7266141051993671, 20.4274178936896966)
#  above was with noIC=TRUE
  tst <-  c(12.2212706117814953,  21.2410851409922579)

 if( fuzz < max(abs( mc100SummaryML$meanSDhatf  - tst ))) {
    cat("Calculated value is not the same as test value in test 21. Value:\n")
    printTestValue(mc100SummaryML$meanSDhatf, digits=18)
    cat("difference:\n")
    print(mc100SummaryML$meanSDhatf - tst, digits=18)
    all.ok <- FALSE  
    }

# using early version of GPA 
# tst <- t(matrix(c( 
#      -0.43072700321218, -0.669882007666349,
#      -2.45915806739033,  0.443541199032001,
#      -0.39598579113561, -0.055347774362779,
#      -1.28804986891355, -2.462726422671402,
#       1.03579550708264,  0.413078009465796,
#       2.50999408754555, -0.359689824353282), 2,6))

 tst <- t(matrix(c(
      -0.430727502635209092, -0.669890058809444078 , 
      -2.45916993900246084 ,  0.443526818436103198 , 
      -0.395990134050530962, -0.0553539782496834665, 
      -1.28804570353536718 , -2.46275799746584845  , 
       1.0358428846874721  , 0.413067269627077849  , 
       2.51004726054249838 , -0.359679178084249429), 2,6))
#  above was with noIC=TRUE
 tst <- t(matrix(c(
      -0.632879784699293069 , -0.354831582672403556 , 
      -1.59531543933821141 , 0.582831484242154474 , 
      -0.399060012837218814 , -0.0892344631765265017 , 
      -1.83055582409809148 , -0.331659929649745777 , 
      1.15013393059134073 , -1.14983030160259503 , 
      0.402743877294618802 , 1.06826205318063217 ), 2,6)) # Linux Mandrake 32bit
 tst <- t(matrix(c( 
      -0.632879702135179301 , -0.354831965108128422 , 
      -1.59531627367130113 , 0.582832449049474732 , 
      -0.399060005892669878 , -0.0892348968104166307 , 
      -1.83055540872476286 , -0.331661458302995982 , 
      1.15013525749168988 , -1.14983284599435365 , 
      0.402743928085540492 , 1.06826581647906949), 2,6))# Solaris
fuzz <- 1e-5  

 if( fuzz < max(abs( mc100SummaryML$meanhatB.error - tst ))) {
    cat("Calculated value is not the same as test value in test 23. Value:\n")
    printTestValue(mc100SummaryML$meanhatB.error, digits=18)
    cat("difference:\n")
    print(mc100SummaryML$meanhatB.error - tst, digits=18)
    all.ok <- FALSE  
    }

# using early version of GPA 
# tst <- t(matrix(c( 
#	 2.1177074120071,  2.8630417742169,
#	 5.8729782096725,  8.9365319798381,
#	 1.4920603854630,  1.9177841472202,
#	 6.9114403738557, 11.0123879149370,
#	 9.6088971180058, 10.4635264914877,
#	12.4426965380737, 17.4104141305891), 2,6))

 tst <- t(matrix(c( 
      2.11774533901723538,  2.86303620239022649, 
      5.87290229669339681,  8.9365652774522939 , 
      1.49204785812003604,  1.91779395295759936, 
      6.91158798109280426, 11.0123561946326376 , 
      9.60890954609334536, 10.463541675796586  , 
     12.4426951488466564,  17.4103921927887733 ),2,6))
#  above was with noIC=TRUE
 tst <- t(matrix(c( 
      2.33840665475836795 , 3.13852181903253502 , 
      6.66951161084711774 , 8.71483640509499935 , 
      1.33929845112549506 ,  2.0031133449915548 , 
      8.52324682951719304 , 11.2751832187060881 , 
      10.5144095556300563 , 7.63423805298109581 , 
      12.7681497781142195 , 15.1413132072165837 ),2,6))# Linux Mandrake 32bit
 tst <- t(matrix(c( 
      2.33840667415145909 , 3.13852213524006274 , 
      6.66951154478138974 ,  8.7148370261171344 , 
      1.33929845221169086 , 2.00311359935738809 , 
      8.52324720810181091 , 11.2751832991476277 , 
      10.5144106848696381 , 7.63424186177272812 , 
      12.7681497223306177 ,   15.14131310687468),2,6))# Solaris


 if( fuzz < max(abs( mc100SummaryML$SDhatB - tst ))) {
    cat("Calculated value is not the same as test value in test 24. Value:\n")
    printTestValue(mc100SummaryML$SDhatB, digits=18)
    cat("difference:\n")
    print(mc100SummaryML$SDhatB - tst, digits=18)
    all.ok <- FALSE  
    }


 fuzz <- 1e-5

# tst <-  c(0.01709371,        1.211911)
# tst <-  c(0.017095808203623, 1.1266841304661)
# above was when Df was percent change (now PCf) rather than difference

 tst <- c(0.0079738515556172, -0.043994019995693)
#  above was with noIC=TRUE
 tst <- c( -0.00505274415462879378,  -0.042964869565685973)

 if( fuzz < max(abs(  (mc100SummaryML$meanhatDf.error)  - tst ))) {
    cat("Calculated value is not the same as test value in test 27. Value:\n")
    printTestValue(mc100SummaryML$meanhatDf.error, digits=18)
    cat("difference:\n")
    print(mc100SummaryML$meanhatDf.error - tst, digits=18)
    all.ok <- FALSE  
    }

#  tst <- c(1.367774,        94.06914)
#  tst <- c(1.3677449004785, 92.129604373376)
# above was when Df was percent change (now PCf) rather than difference
 tst <- c( 0.67484628813627, 1.0077586888048)
#  above was with noIC=TRUE
 tst <- c(0.696639940140005298,  0.994209984498502064)

 if( fuzz < max(abs( mc100SummaryML$meanSDhatDf  - tst ))) {
    cat("Calculated value is not the same as test value in test 29. Value:\n")
    printTestValue(mc100SummaryML$meanSDhatDf, digits=18)
    cat("difference:\n")
    print(mc100SummaryML$meanSDhatDf - tst, digits=18)
    all.ok <- FALSE  
    }

cat("tests completed.\n")

if (! all.ok) stop("some tests FAILED")
} else cat("CDNmoney data not available. Tests skipped.\n")
