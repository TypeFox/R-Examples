.tegarchRecursion <- function (iN, omega, phi1, kappa1, kappastar, df, skew2, dfpluss1, 
    mueps, y, signnegy, signarg, skewterm, lambda, lambdadagg, 
    u) 
.C("tegarchrecursion", iN = as.integer(iN), omega = as.double(omega), phi1 = as.double(phi1), 
    kappa1 = as.double(kappa1), kappastar = as.double(kappastar), 
    df = as.double(df), skew2 = as.double(skew2), dfpluss1 = as.double(dfpluss1), 
    mueps = as.double(mueps), y = as.double(y), signnegy = as.double(signnegy), 
    signarg = as.double(signarg), skewterm = as.double(skewterm), 
    lambda = as.double(lambda), lambdadagg = as.double(lambdadagg), 
    u = as.double(u), PACKAGE = "betategarch")

.tegarchRecursion2 <- function (iN, omega, phi1, phi2, kappa1, kappa2, kappastar, df, 
    skew2, dfpluss1, mueps, y, y2, signnegy, signarg, skewterm, 
    lambda, lambda1dagg, lambda2dagg, u) 
.C("tegarchrecursion2", iN = as.integer(iN), omega = as.double(omega), phi1 = as.double(phi1), 
    phi2 = as.double(phi2), kappa1 = as.double(kappa1), kappa2 = as.double(kappa2), 
    kappastar = as.double(kappastar), df = as.double(df), skew2 = as.double(skew2), 
    dfpluss1 = as.double(dfpluss1), mueps = as.double(mueps), 
    y = as.double(y), y2 = as.double(y2), signnegy = as.double(signnegy), 
    signarg = as.double(signarg), skewterm = as.double(skewterm), 
    lambda = as.double(lambda), lambda1dagg = as.double(lambda1dagg), 
    lambda2dagg = as.double(lambda2dagg), u = as.double(u), PACKAGE = "betategarch")
