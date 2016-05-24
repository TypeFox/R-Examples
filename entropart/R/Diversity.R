Diversity <-
function (NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("Diversity")
}


Diversity.ProbaVector <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
  
  Entropy <- Tsallis(NorP, q=q, CheckArguments=FALSE)
  
  return (expq(Entropy, q))
}


Diversity.AbdVector <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return (bcDiversity(Ns=NorP, q=q, Correction=Correction, CheckArguments=CheckArguments))
}


Diversity.integer <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP Ns must be provided.")
    }
  }
  return (bcDiversity(Ns=NorP, q=q, Correction=Correction, CheckArguments=CheckArguments))
}


Diversity.numeric <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
    return (Diversity.ProbaVector(NorP, q=q, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return (Diversity.AbdVector(NorP, q=q, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcDiversity <-
function(Ns, q = 1, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  Entropy <- bcTsallis(Ns, q, Correction, CheckArguments=FALSE)
  
  return (expq(Entropy, q))
}
