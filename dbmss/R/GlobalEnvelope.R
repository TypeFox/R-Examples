GlobalEnvelope <-
function(Simulations, Alpha) {
  
  verifyclass(Simulations, "fv")  
  
  # Initialize
  Simulations <- as.data.frame(Simulations)[, -1] # Eliminate r
  NumberOfSimulations <- ncol(Simulations)
  TargetNumber <- (1-Alpha)*NumberOfSimulations
  KeptSimulations <- Simulations
  PresentMinValues <- apply(Simulations, 1, min, na.rm=TRUE)
  PresentMaxValues <- apply(Simulations, 1, max, na.rm=TRUE)
  # Loop until the target number of simulations is kept
  while(ncol(KeptSimulations) > TargetNumber) {
    # Remember previous min and max
    PreviousMinValues <- PresentMinValues
    PreviousMaxValues <- PresentMaxValues
    # Select the simulations that gave extreme values
    SimulationsToDrop <- c(unlist(apply(KeptSimulations, 1, which.min)), unlist(apply(KeptSimulations, 1, which.max)))
    # Drop them
    KeptSimulations <- KeptSimulations[, -SimulationsToDrop]  
    # Fails if no simulations are left
    if (is.null(dim(KeptSimulations)))
      stop("Global envelope could not be calculated. More simulations are necessary.")
    # Calculate min and max
    PresentMinValues <- apply(KeptSimulations, 1, min, na.rm=TRUE)
    PresentMaxValues <- apply(KeptSimulations, 1, max, na.rm=TRUE)
  }  
  # Interpolate because the kept number of simulations is not always the target
  NumberOfKeptSimulations <- ncol(KeptSimulations)
  Glo <- PresentMinValues + (PreviousMinValues-PresentMinValues)/NumberOfSimulations*(TargetNumber-NumberOfKeptSimulations)
  Ghi <- PresentMaxValues + (PreviousMaxValues-PresentMaxValues)/NumberOfSimulations*(TargetNumber-NumberOfKeptSimulations)
  return(rbind(Glo, Ghi))
}
