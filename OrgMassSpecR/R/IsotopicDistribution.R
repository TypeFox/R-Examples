IsotopicDistribution <- function(formula = list(), charge = 1) {

# The formula should include the hydrogens added or subracted by the charge.
  
  # Error if charge = 0
  if(charge == 0) stop("a charge of zero is not allowed")

  inputFormula <- list(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0, Br = 0, Cl = 0, F = 0, Si = 0)
  inputFormula[names(formula)] <- formula  # replace default values with argument values
  
  simulation <- function(inputFormula) {
    
    # Simulate isotopic distribution of one molecule.
    
    massCarbon <- sum(sample(c(12.0000000, 13.0033548378),
      size = inputFormula$C,
      replace = TRUE,
      prob = c(0.9893, 0.0107)))
    
    massHydrogen <- sum(sample(c(1.0078250321, 2.0141017780),
      size = inputFormula$H,
      replace = TRUE,
      prob = c(0.999885, 0.000115)))
     
    massNitrogen <- sum(sample(c(14.0030740052, 15.0001088984),
      size = inputFormula$N,
      replace = TRUE,
      prob = c(0.99632, 0.00368)))
      
    massOxygen <- sum(sample(c(15.9949146221, 16.99913150, 17.9991604),
      size = inputFormula$O,
      replace = TRUE,
      prob = c(0.99757, 0.00038, 0.00205)))
    
    massSulfer <- sum(sample(c(31.97207069, 32.97145850, 33.96786683, 35.96708088),
      size = inputFormula$S,
      replace = TRUE,
      prob = c(0.9493, 0.0076, 0.0429, 0.0002)))
    
    massPhosphorus <- inputFormula$P * 30.97376151
    
    massBromine <- sum(sample(c(78.9183376, 80.916291),
      size = inputFormula$Br,
      replace = TRUE,
      prob = c(0.5069, 0.4931)))
    
    massChlorine <- sum(sample(c(34.96885271, 36.96590260),
      size = inputFormula$Cl,
      replace = TRUE,
      prob = c(0.7578, 0.2422)))
    
    massFluorine <- inputFormula$F * 18.99840320
    
    massSilicon <- sum(sample(c(27.9769265327, 28.97649472, 29.97377022),
      size = inputFormula$Si,
      replace = TRUE,
      prob = c(0.922297, 0.046832, 0.030872)))
    
    massMolecule <- sum(massCarbon, massHydrogen, massNitrogen, massOxygen,
                        massSulfer, massPhosphorus, massBromine, massChlorine,
                        massFluorine, massSilicon)
    
    mz <- massMolecule / abs(charge)  # m/z of molecule       
              
    return(mz)
  }
  
  # Simulate set of molecules.
  sim <- replicate(10000, expr = simulation(inputFormula))
  
  # Bin molecules according to their m/z value and make spectrum.
  b <- seq(from = min(sim) - (1 / (2 * abs(charge))),
            to = max(sim) + 1, 
            by = 1 / abs(charge))
  bins <- cut(sim, breaks = b)
  mz <- round(tapply(sim, bins, mean), digits = 2)
  intensity <- as.vector(table(bins))
  spectrum <- data.frame(mz, intensity)
  spectrum <- spectrum[spectrum$intensity != 0, ]
  spectrum$percent <- with(spectrum, round(intensity / max(intensity) * 100, digits = 2))
  row.names(spectrum) <- 1:(nrow(spectrum))
      
  return(spectrum)

}
