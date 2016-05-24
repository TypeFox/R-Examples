MonoisotopicMass <- function(formula = list(), isotopes = list(), charge = 0) {
    
    defaultFormula <- list(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0, Br = 0, Cl = 0, F = 0, Si = 0, x = 0)
    defaultFormula[names(formula)] <- formula   # replace default values with argument values
    
    defaultIsotopes <- list(C = 12, 
                            H = 1.0078250321, 
                            N = 14.0030740052, 
                            O = 15.9949146221, 
                            S = 31.97207069, 
                            P = 30.97376151,
                            Br = 78.9183376,
                            Cl = 34.96885271,
                            F = 18.99840320,
                            Si = 27.9769265327,
                            x = 0)

    defaultIsotopes[names(isotopes)] <- isotopes

    if(charge < 0 & abs(charge) > defaultFormula$H)
        stop("the number of negative charges exceeds the number of hydrogens in the formula list")

    mass <- (defaultFormula$C * defaultIsotopes$C + 
             defaultFormula$H * defaultIsotopes$H +
             defaultFormula$N * defaultIsotopes$N + 
             defaultFormula$O * defaultIsotopes$O +
             defaultFormula$S * defaultIsotopes$S + 
             defaultFormula$P * defaultIsotopes$P +
             defaultFormula$Br * defaultIsotopes$Br +
             defaultFormula$Cl * defaultIsotopes$Cl +
             defaultFormula$F * defaultIsotopes$F +
             defaultFormula$Si * defaultIsotopes$Si +
             defaultFormula$x * defaultIsotopes$x)

    if(charge != 0) mass <- abs((mass + charge * 1.007276466) / charge)

    return(mass)

}
