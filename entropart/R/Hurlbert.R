Hurlbert <-
function(NorP, k = 2, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("Hurlbert")
}


Hurlbert.ProbaVector <-
function(NorP, k = 2, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
  
  index <- sum(1 - (1-NorP)^k)
  names(index) <- "Biased"
  return (index)
}

Hurlbert.AbdVector <-
function(NorP, k = 2, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return (bcHurlbert(Ns=NorP, k=k, CheckArguments=CheckArguments))
}


Hurlbert.integer <-
function(NorP, k = 2, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return (bcHurlbert(Ns=NorP, k=k, CheckArguments=CheckArguments))
}


Hurlbert.numeric <-
function(NorP, k = 2, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
    return (Hurlbert.ProbaVector(NorP, k=k, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return (Hurlbert.AbdVector(NorP, k=k, CheckArguments=CheckArguments))
  }
}


bcHurlbert <-
function(Ns, k = 2, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  N <- sum(Ns)
  if (k > N) stop("The order of diversity cannot be greater than the size of the sample (check that k <= sum(Ns)).")
  S <- length(Ns)
  # Use lchoose and differences to avoid Inf
  lcNk <- lchoose(N, k)
  index <- S - sum(exp(lchoose(N-Ns, k)-lcNk))
  names(index) <- "Unbiased"
  return (index)
  
}
