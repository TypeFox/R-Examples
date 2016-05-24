Rao <-
function(NorP, Tree, Correction = "Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("Rao")
}


Rao.ProbaVector <-
function(NorP, Tree, Correction = "Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
  
  return (AllenH(NorP, q=2, PhyloTree=Tree, Normalize=FALSE, CheckArguments=FALSE))
}


Rao.AbdVector <-
function(NorP, Tree, Correction = "Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return (bcRao(Ns=NorP, Tree=Tree, Correction=Correction, CheckArguments=CheckArguments))
}


Rao.integer <-
function(NorP, Tree, Correction = "Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return (bcRao(Ns=NorP, Tree=Tree, Correction=Correction, CheckArguments=CheckArguments))
}


Rao.numeric <-
function(NorP, Tree, Correction = "Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
    return (Rao.ProbaVector(NorP, Tree=Tree, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return (Rao.AbdVector(NorP, Tree=Tree, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcRao <-
function(Ns, Tree, Correction="Lande", CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()

  # Eliminate 0
  Ns <- Ns[Ns > 0]
  N <- sum(Ns)
  # Exit if Ns contains no or a single species
  if (length(Ns) < 2) {
  	if (length(Ns) == 0) {
  	  entropy <- NA
  	  names(entropy) <- "No Species"
  	  return (entropy)
  	} else {
  	  entropy <- 0
  	  names(entropy) <- "Single Species"
  	  return (entropy)
  	}
  } else {
    # Probabilities instead of abundances
    if (N < 2) {
      warning("Bias correction attempted with probability data. Correction forced to 'None'")
      Correction <- "None"
    }
  }	
 
  if (Correction == "Lande"  | Correction == "Best") {
    N <- sum(Ns)
    entropy <- N/(N-1)*AllenH(Ns/sum(Ns), q=2, PhyloTree=Tree, Normalize=FALSE, CheckArguments=FALSE)
    names(entropy) <- Correction
    return (entropy) 
  } else {
    return (bcPhyloEntropy(Ns, q = 2, Tree = Tree, Normalize = FALSE, Correction = Correction, CheckArguments = FALSE)$Total) 
  }
}
