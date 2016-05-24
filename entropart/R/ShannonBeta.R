ShannonBeta <-
function(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  UseMethod("ShannonBeta")
}


ShannonBeta.ProbaVector <-
function(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (TsallisBeta(NorP, NorPexp, q=1, CheckArguments=FALSE))
}


ShannonBeta.AbdVector <-
function(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
{
  return (bcShannonBeta(Ns=NorP, Nexp=NorPexp, Correction=Correction, CheckArguments=CheckArguments))
}


ShannonBeta.integer <-
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
  return (bcShannonBeta(Ns=NorP, Nexp=NorPexp, Correction=Correction, CheckArguments=CheckArguments))
}


ShannonBeta.numeric <-
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
    return (ShannonBeta.ProbaVector(NorP, NorPexp, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return (ShannonBeta.AbdVector(NorP, NorPexp, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcShannonBeta <-
function(Ns, Nexp, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (Correction == "ZhangGrabchak") {
    # Eliminate 0
    Y <- Nexp[Ns > 0]
    X <- Ns[Ns > 0]
    m <- sum(Nexp)
    n <- sum(Ns)
    Ps <- X/n
    V1 <- 1:m
    V2 <- 1:(n-1)
    # p_V_Ns1 is an array, containing (1 - Y_s/(m-j+1)) for each species (lines) and all j from 1 to m (Y_s may be 0)
    p_V_Ns1 <- outer(Y, V1, function(Y, j) 1- Y/(m-j+1))
    # p_V_Ns2 is an array, containing (1 - (X_s-1)/(n-j)) for each species (lines) and all j from 1 to n-1
    p_V_Ns2 <- outer(X, V2, function(X, j) 1- (X-1)/(n-j))
    # Useful values are products from j=1 to v, so prepare cumulative products
    p_V_Ns1 <- t(apply(p_V_Ns1, 1, cumprod))
    p_V_Ns2 <- t(apply(p_V_Ns2, 1, cumprod))
    # Sum of products weighted by w_v=1/v
    w_v1 <- 1/V1
    w_v2 <- 1/V2
    S_v <- function(s) {
      Usedv1 <- 1:(m-Y[s])
      Usedv2 <- 1:(n-X[s])
      return (sum(w_v1[Usedv1]*p_V_Ns1[s, Usedv1]) - sum(w_v2[Usedv2]*p_V_Ns2[s, Usedv2]))
    }
    entropy <- sum(Ps*vapply(1:length(Ps), S_v, 0))
    names(entropy) <- Correction
    return (entropy)
  } else {
    return (bcTsallisBeta(Ns, Nexp, 1, Correction, CheckArguments=FALSE))    
  }
}
