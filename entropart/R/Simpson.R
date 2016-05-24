Simpson <-
function(NorP, Correction="Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("Simpson")
}


Simpson.ProbaVector <-
function(NorP, Correction="Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ps)) {
      NorP <- Ps
    } else {
      stop("An argument NorP or Ps must be provided.")
    }
  }
  if (CheckArguments)
    CheckentropartArguments()
  
  return (Tsallis.ProbaVector (NorP, q=2, CheckArguments=FALSE))
}


Simpson.AbdVector <-
function(NorP, Correction="Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return (bcSimpson(Ns=NorP, Correction=Correction, CheckArguments=CheckArguments))
}


Simpson.integer <-
function(NorP, Correction="Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return (bcSimpson(Ns=NorP, Correction=Correction, CheckArguments=CheckArguments))
}


Simpson.numeric <-
function(NorP, Correction="Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ps)) {
      NorP <- Ps
    } else {
      if (!missing(Ns)) {
        NorP <- Ns
      } else {
        stop("An argument NorP or Ps or Ns must be provided.")
      }
    }
  }
  
  if (abs(sum(NorP) - 1) < length(NorP)*.Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error
    return (Simpson.ProbaVector(NorP, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return (Simpson.AbdVector(NorP, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcSimpson <-
function(Ns, Correction="Lande", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()

  # Eliminate 0
  Ns <- Ns[Ns > 0]
  N <- sum(Ns)
  # Exit if Ns contains no or a single species
  if (length(Ns) < 2) {
  	if (length(Ns) == 0) {
  		return(NA)
  	} else {
  		return(0)
  	}
  } else {
    # Probabilities instead of abundances
    if (N < 2) {
      warning("Bias correction attempted with probability data. Correction forced to 'None'")
      Correction <- "None"
    }
  }

  if (Correction == "Lande" | Correction == "Best") {
    entropy <- N/(N-1)*Tsallis(as.ProbaVector(Ns), 2, CheckArguments=FALSE)
    names(entropy) <- Correction
    return (entropy)
  } else {
    return (bcTsallis(Ns, q=2, Correction, CheckArguments=FALSE)) 
  }
}
