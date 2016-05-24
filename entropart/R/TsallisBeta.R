TsallisBeta <-
function(NorP, NorPexp = NULL, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  UseMethod("TsallisBeta")
}


TsallisBeta.ProbaVector <-
function(NorP, NorPexp = NULL, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (length(NorP) != length(NorPexp)) {
    stop("NorP and NorPexp should have the same length.")
  }  
  
  dataBeta <- NorP^q * lnq(NorP/NorPexp, q)
  dataBeta[NorP == 0] <- 0
  
  entropy <- sum(dataBeta)
  names(entropy) <- "None"
  return (entropy)
}


TsallisBeta.AbdVector <-
function(NorP, NorPexp = NULL, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  return (bcTsallisBeta(Ns=NorP, Nexp=NorPexp, q=q, Correction=Correction, CheckArguments=CheckArguments))
}


TsallisBeta.integer <-
function(NorP, NorPexp = NULL, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
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
  return (bcTsallisBeta(Ns=NorP, Nexp=NorPexp, q=q, Correction=Correction, CheckArguments=CheckArguments))
}


TsallisBeta.numeric <-
function(NorP, NorPexp = NULL, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
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
    return (TsallisBeta.ProbaVector(NorP, NorPexp, q=q, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return (TsallisBeta.AbdVector(NorP, NorPexp, q=q, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcTsallisBeta <-
function(Ns, Nexp = NULL, q, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (length(Ns) != length(Nexp)) {
    stop("Ns and Nexp should have the same length.")
  }  
  
  # No correction
  if (Correction == "None") {
    return (TsallisBeta.ProbaVector(Ns/sum(Ns), Nexp/sum(Nexp), q, CheckArguments=FALSE))
  }
  
  # Sample coverage
  Nrecords <- sum(Ns)
  SampleCoverage <- Coverage(Ns, CheckArguments=FALSE)
  # Sample coverage (expected)
  Nrecordsexp <- sum(Nexp)
  SampleCoverageexp <- Coverage(Nexp, CheckArguments=FALSE)
  if (Correction == "ChaoShen" | Correction == "Best") {
    CiPsi <- SampleCoverage * Ns / Nrecords
    CPs <- SampleCoverageexp * Nexp / Nrecordsexp
    dataBeta <- CiPsi^q * lnq(CiPsi/CPs, q) / (1 -(1-CiPsi)^Nrecords)    
    # force 0log0=0                                                                                                                                       
    dataBeta[Ns == 0] <- 0
    entropy <- sum(dataBeta)
    names(entropy) <- Correction
    return (entropy)
  } 
  warning("Correction was not recognized")
  return (NA)
}
