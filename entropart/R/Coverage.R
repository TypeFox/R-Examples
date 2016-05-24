Coverage <-
function(Ns, Estimator = "ZhangHuang", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Round values
  Ns <- as.integer(round(Ns))
  # Eliminate zeros
  Ns <- Ns[Ns>0]
  # Calculate abundance distribution
  DistN <- tapply(Ns, Ns, length)
  Singletons <- DistN["1"]

  # No singletons, C=1
  if (is.na(Singletons)) {
    coverage <- 1
    names(coverage) <- "No singleton"
    return (coverage)
  }
  
  # Abundances
  Nu <- as.integer(names(DistN))
  SampleSize <- sum(Ns)

  # Singletons only
  if (Singletons == SampleSize) {
    warning ("Sample coverage is 0, most bias corrections will return NaN.")
    coverage <- 0
    names(coverage) <- "Singletons only"
    return (coverage)
  }

  if (Estimator == "ZhangHuang") {
    Ps <- Ns/SampleSize
    if (any(Ps >= .5)) {
      warning ("Zhang-Huang sample coverage cannot be estimated because one probability is over 1/2. Chao estimator is returned.")
      Estimator <- "Chao"
    } else {
      # Use Nu%%2*2-1 for (-1)^(Nu+1)
      coverage <- 1 - sum((Nu%%2*2-1) / choose(SampleSize, Nu) * DistN)
      names(coverage) <- Estimator
      return (coverage)
    }    
  }
  if (Estimator == "Chao") {
    Doubletons <- DistN["2"]
    if (is.na(Doubletons)) {
      warning ("Chao's sample coverage cannot be estimated because there are no doubletons. Turing estimator is returned.")
      Estimator <- "Turing"
    } else {
      coverage <- 1 - Singletons / SampleSize *((SampleSize - 1) * Singletons / ((SampleSize - 1) * Singletons + 2 * Doubletons))
      names(coverage) <- Estimator
      return (coverage)
    }
  }
  if (Estimator == "Turing") {
    coverage <- 1 - Singletons / SampleSize
    names(coverage) <- Estimator
    return (coverage)
  }
  warning("Correction has not been recognized")
  return (NA)

}
