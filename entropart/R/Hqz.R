Hqz <-
function(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("Hqz")
}


Hqz.ProbaVector <-
function(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
  
  Ps <- NorP
  # If names are missing, the probability vector and the similarity vector are assumed to be in the same order
  if (is.null(colnames(Z)) | is.null(names(Ps))) {
    if (ncol(as.matrix(Z)) != length(Ps)) # as.matrix(Z) in case it has been reduced to a single value because of zeros
      # The matrix is square (this has been checked)
      stop("The matrix dimension must equal the probability vector length.")    
    # Eliminate zeros
    Z <- as.matrix(Z)[Ps != 0, Ps != 0]
    Ps <- Ps[Ps != 0]
  } else { # Matrix and Ps are be named.  
    # Eliminate zeros
    Ps <- Ps[Ps != 0]
    if (length(setdiff(names(Ps), colnames(Z))) != 0)
      # The matrix is square (this has been checked)
      stop("Some species are missing in the similarity matrix.")    
    Z <- as.matrix(Z)[names(Ps), names(Ps)]
  }
  
  # Calculate (Zp)
  Zp <- Z %*% Ps
  if (q == 1) {
    # Limit value
    Entropy <- -Ps %*% log(Zp)
  } else {
    # Calculate (Zp)^(q-1)
    Zpqm1 <- Zp^(q-1)
    # Calculate Hqz
    Entropy <- (1-(Ps %*% Zpqm1))/(q-1)
  }
  # Return the value of entropy, as a number rather than a 1x1 matrix
  Entropy <- as.numeric(Entropy)
  names(Entropy) <- "None"
  return (Entropy)
}


Hqz.AbdVector <-
function(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcHqz(Ns=NorP, q=q, Z=Z, Correction=Correction, CheckArguments=CheckArguments))
}


Hqz.integer <-
function(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcHqz(Ns=NorP, q=q, Z=Z, Correction=Correction, CheckArguments=CheckArguments))
}


Hqz.numeric <-
function(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
    return(Hqz.ProbaVector(NorP, q=q, Z=Z, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return(Hqz.AbdVector(NorP, q=q, Z=Z, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcHqz <- 
function (Ns, q = 1, Z = diag(length(Ns)), Correction = "Best", CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  # If names are missing, the probability vector and the similarity vector are assumed to be in the same order
  if (is.null(colnames(Z)) | is.null(names(Ns))) {
    if (ncol(as.matrix(Z)) != length(Ns)) # as.matrix(Z) in case it has been reduced to a single value because of zeros
      # The matrix is square (this has been checked)
      stop("The matrix dimension must equal the abundance vector length.")    
    # Eliminate zeros
    Z <- as.matrix(Z)[Ns != 0, Ns != 0]
    Ns <- Ns[Ns != 0]
  } else { # Matrix and Ns are be named.  
    # Eliminate zeros
    Ns <- Ns[Ns != 0]
    if (length(setdiff(names(Ns), colnames(Z))) != 0)
      # The matrix is square (this has been checked)
      stop("Some species are missing in the similarity matrix.")    
    Z <- as.matrix(Z)[names(Ns), names(Ns)]
  }
  
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
  Ps <- Ns/N

  
  # No correction
  if (Correction == "None") {
    return (Hqz.ProbaVector(Ps, q, Z, CheckArguments = FALSE))
  }
  
  if (Correction == "MarconZhang" | Correction == "Best") {
    V <- 1:(N-1)
    # p_V_Ns is an array, containing (1 - (n_s-1)/(n-j)) for each species (lines) and all j from 1 to n-1
    p_V_Ns <- outer(Ns, V, function(Ns, j) 1- (Ns-1)/(N-j))
    # Useful values are products from j=1 to v, so prepare cumulative products
    p_V_Ns <- apply(p_V_Ns, 1, cumprod)
    # Sum of products weighted by w_v
    S_v <- function(s) {
      Usedv <- 1:(N-Ns[s])
      return (sum(w_v[Usedv]*p_V_Ns[Usedv, s]))
    }
  }
  
  # Sample coverage
  C <- Coverage(Ns)
  CPs <- C*Ps
  # Calculate the average similarity
  Zcopy <- Z
  diag(Zcopy) <- NA
  AverageZ <- mean(Zcopy, na.rm=TRUE)
  # Add (1-C)*AverageZ to Zp
  Zp <- as.vector(Z %*% CPs) + rep(AverageZ*(1-C), length(CPs))
  
  if (Correction == "ChaoShen" | Correction == "Best") {
    # Horvitz-Thomson multiplier (replaces Ps)
    HTCPs <- CPs/(1 - (1-CPs)^N)
    # Force 0/0=0 and 0log0=0
    HTCPs[CPs == 0] <- 0
    Zp[Zp == 0] <- 1
    HT <- (sum(HTCPs*lnq(1/Zp, q)))
  }
  
  if ((Correction == "MarconZhang" | Correction == "Best") & (q != 1)) {
    Zpqm1 <- Zp^(q-1)
    # Force 0^(q-1)=0
    Zpqm1[Zp == 0] <- 0
    K <- sum(CPs * Zpqm1)
    # Weights
    i <- 1:N
    w_vi <- (1-AverageZ)*(i-q)/i
    w_v <- cumprod(w_vi)
    Taylor <- 1 + sum(Ps*vapply(1:length(Ns), S_v, 0))
    FirstTerms <- CPs*(AverageZ+(1-AverageZ)*CPs)^(q-1)
    U <- Taylor-sum(FirstTerms)
    MZ <- ((K+U-1)/(1-q))
  }
  
  if ((Correction == "MarconZhang" | Correction == "Best") & (q == 1)) {
    L <- -sum(CPs*log(Zp))
    # Weights
    w_v <- ((1-AverageZ)^V)/V
    Taylor <- sum(Ps*vapply(1:length(Ns), S_v, 0))
    FirstTerms <- -CPs*log(AverageZ+(1-AverageZ)*CPs)
    X <- Taylor-sum(FirstTerms)
    # MZ <- -sum(CPs*log(Zp)) -(1-C)*log(AverageZ)
    MZ <- L+X
  }
  
  
  if (Correction == "ChaoShen") {
    names(HT) <- Correction
    return(HT)
  }
  if (Correction == "MarconZhang") {
    names(MZ) <- Correction
    return(MZ)
  }
  if (Correction == "Best") {
    entropy <- max(HT, MZ)
    names(entropy) <- Correction
    return(entropy)
  }
  warning("Correction was not recognized")
  return (NA)
}
