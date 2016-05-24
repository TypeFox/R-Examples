GenSimpson <-
  function(NorP, r = 1, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
  {
    UseMethod("GenSimpson")
  }


GenSimpson.ProbaVector <-
function(NorP, r = 1, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
  
  entropy <- sum(sum(NorP*(1-NorP)^r))
  names(entropy) <- "Biased"
  return (entropy)
}


GenSimpson.AbdVector <-
function(NorP, r = 1, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return (bcGenSimpson(Ns=NorP, r=r, CheckArguments=CheckArguments))
}


GenSimpson.integer <-
  function(NorP, r = 1, CheckArguments = TRUE, Ps = NULL, Ns = NULL)
  {
    if (missing(NorP)){
      if (!missing(Ns)) {
        NorP <- Ns
      } else {
        stop("An argument NorP or Ns must be provided.")
      }
    }
    return(bcGenSimpson(Ns=NorP, r=r, CheckArguments=CheckArguments))
  }


GenSimpson.numeric <-
  function(NorP, r = 1, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
      return(GenSimpson.ProbaVector(NorP, r=r, CheckArguments=CheckArguments))
    } else {
      # Abundances
      return(GenSimpson.AbdVector(NorP, r=r, CheckArguments=CheckArguments))
    }
  }


bcGenSimpson <-
function(Ns, r = 1, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  
  entropy <- EntropyEstimation::GenSimp.z(Ns, r)
  names(entropy) <- "Unbiased"
  return (entropy)
}
