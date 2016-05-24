HqzBeta <-
function(NorP, NorPexp = NULL, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  UseMethod("HqzBeta")
}


HqzBeta.ProbaVector <-
function(NorP, NorPexp = NULL, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  Ps <- NorP
  Pexp <- NorPexp
  
  if (length(Ps) != length(Pexp)) {
    stop("Ps and Pexp should have the same length.")
  }  
  # If names are missing, the probability vector and the similarity vector are assumed to be in the same order
  if (is.null(colnames(Z)) | is.null(names(Ps))) {
    if (ncol(as.matrix(Z)) != length(Ps))  # as.matrix(Z) in case it has been reduced to a single value because of zeros
      # The matrix is square (this has been checked)
      stop("The matrix dimension must equal the probability vector length.")    
    # Eliminate zeros
    Z <- as.matrix(Z)[Ps != 0, Ps != 0]
    Pexp <- Pexp[Ps != 0]
    Ps <- Ps[Ps != 0]
  } else { # Matrix and Ps are be named.
    # Reorder Ps and Pexp
    if (!setequal(names(Ps), names(Pexp)))
      stop("Ps and Pexp should have the names.")
    Pexp <- Pexp[names(Ps)]
    # Eliminate zeros
    Pexp <- Pexp[Ps != 0]
    Ps <- Ps[Ps != 0]
    if (length(setdiff(names(Ps), colnames(Z))) != 0)
      # The matrix is square (this has been checked)
      stop("Some species are missing in the similarity matrix.")    
    Z <- as.matrix(Z)[names(Ps), names(Ps)]
  }
  
  # Calculate (Zp)
  Zps <- Z %*% Ps
  Zpexp <- Z %*% Pexp
  
  dataBeta <- Ps * (lnq(1/Zpexp, q)-lnq(1/Zps, q))
  entropy <- sum(dataBeta)
  names(entropy) <- "None"
  return (entropy)
}


HqzBeta.AbdVector <-
function(NorP, NorPexp = NULL, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  return (bcHqzBeta(Ns=NorP, Nexp=NorPexp , q=q, Z=Z, Correction=Correction, CheckArguments=CheckArguments))
}


HqzBeta.integer <-
function(NorP, NorPexp = NULL, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
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
  return (bcHqzBeta(Ns=NorP, Nexp=NorPexp, q=q, Z=Z, Correction=Correction, CheckArguments=CheckArguments))
}


HqzBeta.numeric <-
function(NorP, NorPexp = NULL, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
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
    return (HqzBeta.ProbaVector(NorP, NorPexp, q=q, Z=Z, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return (HqzBeta.AbdVector(NorP, NorPexp, q=q, Z=Z, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcHqzBeta <-
function(Ns, Nexp = NULL, q = 1, Z = diag(length(Ns)), Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (length(Ns) != length(Nexp)) {
    stop("Ns and Nexp should have the same length.")
  }  
  
  # No correction available yet
  if (Correction == "None" | Correction == "Best") {
    return (HqzBeta(Ns/sum(Ns), Nexp/sum(Nexp), q, Z, CheckArguments=FALSE))
  }
  
  warning("Correction was not recognized")
  return (NA)
}
