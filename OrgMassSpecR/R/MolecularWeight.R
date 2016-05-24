MolecularWeight <- function(formula = list(), amu = list()) {
    
    defaultFormula <- list(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0, Br = 0, Cl = 0, F = 0, Si = 0, Sn = 0, x = 0)
    defaultFormula[names(formula)] <- formula   # replace default values with argument values
    
    defaultAmu <- list(C = 12.0107, 
                       H = 1.00794, 
                       N = 14.0067, 
                       O = 15.9994, 
                       S = 32.065, 
                       P = 30.973761, 
                       Br = 79.904,
                       Cl = 35.453,
                       F = 18.9984032,
                       Si = 28.0855,
                       Sn = 118.710,
                       x = 0)

    defaultAmu[names(amu)] <- amu

    mw <- (defaultFormula$C * defaultAmu$C + 
           defaultFormula$H * defaultAmu$H +
           defaultFormula$N * defaultAmu$N + 
           defaultFormula$O * defaultAmu$O +
           defaultFormula$S * defaultAmu$S + 
           defaultFormula$P * defaultAmu$P +
           defaultFormula$Br * defaultAmu$Br + 
           defaultFormula$Cl * defaultAmu$Cl +
           defaultFormula$F * defaultAmu$F +
           defaultFormula$Si * defaultAmu$Si +
           defaultFormula$Sn * defaultAmu$Sn +
           defaultFormula$x * defaultAmu$x)

    return(round(mw, digits = 3))

}
