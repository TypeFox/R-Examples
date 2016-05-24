Dqz <-
function(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("Dqz")
}


Dqz.ProbaVector <-
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
    if (ncol(as.matrix(Z)) != length(Ps))  # as.matrix(Z) in case it has been reduced to a single value because of zeros
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
    Diversity <- exp(-Ps %*% log(Zp))
  } else {
    # Calculate (Zp)^(q-1)
    Zpqm1 <- Zp^(q-1)
    # Calculate Dqz
    Diversity <- (Ps %*% Zpqm1)^(1/(1-q))
  }
  # Return the value of diversity, as a number rather than a 1x1 matrix
  Diversity <- as.numeric(Diversity)
  names(Diversity) <- "None"
  return (Diversity)
}


Dqz.AbdVector <-
function(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcDqz(Ns=NorP, q=q, Z=Z, Correction=Correction, CheckArguments=CheckArguments))
}


Dqz.integer <-
function(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcDqz(Ns=NorP, q=q, Z=Z, Correction=Correction, CheckArguments=CheckArguments))
}


Dqz.numeric <-
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
    return(Dqz.ProbaVector(NorP, q=q, Z=Z, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return(Dqz.AbdVector(NorP, q=q, Z=Z, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcDqz <-
function (Ns, q = 1, Z = diag(length(Ns)), Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  Entropy <- bcHqz(Ns, q, Z, Correction, CheckArguments=FALSE)
  
  return (expq(Entropy, q))
}
