SimpsonBeta <-
function(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  UseMethod("SimpsonBeta")
}


SimpsonBeta.ProbaVector <-
function(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (TsallisBeta(NorP, NorPexp, q=2, CheckArguments=FALSE))
}


SimpsonBeta.AbdVector <-
function(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  return (bcSimpsonBeta(Ns=NorP, Nexp=NorPexp, Correction=Correction, CheckArguments=CheckArguments))
}


SimpsonBeta.integer <-
function(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  if (missing(NorPexp)){
    if (!missing(Nexp)) {
      NorPexp <- Nexp
    } else {
      stop("An argument NorPexp or Nexp must be provided.")
    }
  }
  return (bcSimpsonBeta(Ns=NorP, Nexp=NorPexp, Correction=Correction, CheckArguments=CheckArguments))
}


SimpsonBeta.numeric <-
function(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
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
  if (missing(NorPexp)){
    if (!missing(Pexp)) {
      NorPexp <- Pexp
    } else {
      if (!missing(Nexp)) {
        NorP <- Nexp
      } else {
        stop("An argument NorPexp or Pexp or Nexp must be provided.")
      }
    }
  }
  
  if (abs(sum(NorP) - 1) < length(NorP)*.Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error
    return (SimpsonBeta.ProbaVector(NorP, NorPexp, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return (SimpsonBeta.AbdVector(NorP, NorPexp, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcSimpsonBeta <-
function(Ns, Nexp, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (bcTsallisBeta(Ns, Nexp, 2, Correction, CheckArguments=FALSE))
}
