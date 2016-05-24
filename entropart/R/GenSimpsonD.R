GenSimpsonD <-
function(NorP, r = 1, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("GenSimpsonD")
}


GenSimpsonD.ProbaVector <-
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
  
  diversity <- 1/(1-(sum(NorP*(1-NorP)^r))^(1/r))
  names(diversity) <- "Biased"
  return (diversity)
}


GenSimpsonD.AbdVector <-
function(NorP, r = 1, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return (bcGenSimpsonD(Ns=NorP, r=r, CheckArguments=CheckArguments))
}


GenSimpsonD.integer <-
  function(NorP, r = 1, CheckArguments = TRUE, Ps = NULL, Ns = NULL)
  {
    if (missing(NorP)){
      if (!missing(Ns)) {
        NorP <- Ns
      } else {
        stop("An argument NorP or Ns must be provided.")
      }
    }
    return(bcGenSimpsonD(Ns=NorP, r=r, CheckArguments=CheckArguments))
  }


GenSimpsonD.numeric <-
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
      return(GenSimpsonD.AbdVector(NorP, r=r, CheckArguments=CheckArguments))
    }
  }


bcGenSimpsonD <-
function(Ns, r = 1, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  
  diversity <- 1/(1-(EntropyEstimation::GenSimp.z(Ns, r))^(1/r))
  names(diversity) <- "Unbiased"
  return (diversity)
}
