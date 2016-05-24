Shannon <-
function(NorP, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("Shannon")
}


Shannon.ProbaVector <-
function(NorP, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
  
  return (Tsallis.ProbaVector(NorP, q=1, CheckArguments=FALSE))
}


Shannon.AbdVector <-
function(NorP, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return (bcShannon(Ns=NorP, Correction=Correction, CheckArguments=CheckArguments))
}


Shannon.integer <-
function(NorP, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return (bcShannon(Ns=NorP, Correction=Correction, CheckArguments=CheckArguments))
}


Shannon.numeric <-
function(NorP, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
    return (Shannon.ProbaVector(NorP, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return (Shannon.AbdVector(NorP, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcShannon<-
function(Ns, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (Correction == "Best") Correction <- "ChaoWangJost"
  
  # Eliminate 0
  Ns <- Ns[Ns > 0]
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
  
  # No correction
  if (Correction == "None") {
    return (Shannon.ProbaVector(Ns/sum(Ns), CheckArguments=FALSE))
  }
  
  
  if (Correction == "Miller") {
    return (Shannon(Ns/sum(Ns), CheckArguments=FALSE) + (length(Ns)-1)/2/N)  
  }
  if (Correction == "ChaoShen" | Correction == "GenCov" | Correction == "Marcon") {
    SampleCoverage <- Coverage(Ns, CheckArguments=FALSE)
  }
  if (Correction == "ChaoShen") {
    CPs <- SampleCoverage*Ns/N
  }
  if (Correction == "GenCov" | Correction == "Marcon") {
    CPs <- as.ProbaVector(Ns, Correction="Chao2015", CheckArguments = FALSE)
  } 
  if (Correction == "ChaoShen" | Correction == "GenCov" | Correction == "Marcon") {
    ChaoShen <- sum(-CPs*log(CPs)/(1-(1-CPs)^N))
  }
  if (Correction == "ChaoShen" | Correction == "GenCov") {
    names(ChaoShen) <- Correction
    return (ChaoShen)  
  } 
  if (Correction == "Grassberger") {
    # (-1)^n is problematic for long vectors (returns NA for large values). It is replaced by 1-n%%2*2 (Ns is rounded if is not an integer)
    Grassberger <- sum(Ns/N*(log(N)-digamma(Ns)-(1-round(Ns)%%2*2)/(Ns+1)))
  }
  if (Correction == "Grassberger") {
    names(Grassberger) <- Correction
    return(Grassberger)
  }
  if (Correction == "Marcon") {
    entropy <- max(ChaoShen, Grassberger)
    names(entropy) <- Correction
    return (entropy)
  }
  if (Correction == "Grassberger2003" | Correction == "Schurmann") {
    # Define a function to calculate the integral in the bias formuma for each value of N
    Integral <- function(n, upper) stats::integrate(function(t, n) t^(n-1)/(1+t), 0, upper, n) 
  }
  if (Correction == "Grassberger2003") {
    Integral.V <- unlist(sapply(Ns, Integral, upper = 1)["value",])
  }
  if (Correction == "Schurmann") {
    Integral.V <- unlist(sapply(Ns, Integral, upper = exp(-1/2))["value",])
  }
  if (Correction == "Grassberger2003" | Correction == "Schurmann") {
    entropy <- sum(Ns/N*(digamma(N)-digamma(Ns)-(1-Ns%%2*2)*Integral.V))
    names(entropy) <- Correction
    return (entropy)
  }
  if (Correction == "Holste" | Correction == "Bonachela") {
    seql <- 1:(length(Ns)+N)
    invl <- 1/seql
    cumul <- function(n) {sum(invl[n:length(invl)])}
    suminvl <- vapply(seql, cumul, 0) 
    if (Correction == "Holste") {
      entropy <- sum((Ns+1)/(length(Ns)+N)*suminvl[Ns+2])
      names(entropy) <- Correction
      return (entropy)
    } else {
      entropy <- sum((Ns+1)/(2+N)*suminvl[Ns+2])
      names(entropy) <- Correction
      return (entropy)
    }
  }
  if (Correction == "ChaoWangJost") {
    # Calculate abundance distribution
    DistN <- tapply(Ns, Ns, length)
    Singletons <- DistN["1"]
    if (is.na(Singletons)) Singletons <- 0
    Doubletons <- DistN["2"]
    if (is.na(Doubletons)) Doubletons <- 0
    # Calculate A
    if (Doubletons) {
      A <- 2*Doubletons/((N-1)*Singletons+2*Doubletons)
    } else {
      if (Singletons) {
        A <- 2/((N-1)*(Singletons-1)+2)
      } else {
        A <- 1
      }
    }
    ChaoWangJost <- sum(Ns/N*(digamma(N)-digamma(Ns)))
    if (A != 1) {
      Part2 <- vapply(1:(N-1), function(r) 1/r*(1-A)^r, 0) 
      ChaoWangJost <- as.numeric(ChaoWangJost + Singletons/N*(1-A)^(1-N)*(-log(A)-sum(Part2)))
    }
    names(ChaoWangJost) <- Correction
    return(ChaoWangJost)  
  }
  if (Correction == "ZhangHz") {
    # Values of v
    V <- 1:(N-1)
    Ps <- Ns/N
    # Weight part. Taken in log or goes to Inf for v > 1000; gamma cannot be used for large n, lgamma is preferred.
    lnw_v <- ((V+1)*log(N)+lgamma(N-V)-lgamma(N+1)-log(V))
    # p_V_Ps is an array, containing (1 - p_s - j/n) for each species (lines) and all j from 0 to n-2. Because array indexation starts from 1 in R, j is replaced by j-1.
    p_V_Ps <- outer(Ps, V, function(Ps, j) 1 -Ps -(j-1)/N)
    # Useful values are products from j=0 to v, so prepare cumulative products
    p_V_Ps <- t(apply(p_V_Ps, 1, cumprod))
    # Sum of products, weighted by p_s
    S_s <- function(v) {
      sum(Ps*p_V_Ps[1:length(Ps), v])
    }
    # Apply S_s to all values of v. Use logs or w_v goes to Inf.
    entropy <- sum(exp(lnw_v + log(vapply(V, S_s, 0))))
    names(entropy) <- Correction
    return (entropy)
  }
  if (Correction == "UnveilC") {
    TunedPs <- as.ProbaVector(Ns, Correction="Chao2015", Unveiling="geom", RCorrection = "Chao1", CheckArguments = FALSE)
  }
  if (Correction == "UnveiliC") {
    TunedPs <- as.ProbaVector(Ns, Correction="Chao2015", Unveiling="geom", RCorrection = "iChao1", CheckArguments = FALSE)
  }
  if (Correction == "UnveilJ") {
    TunedPs <- as.ProbaVector(Ns, Correction="Chao2015", Unveiling="geom", RCorrection = "Jackknife", CheckArguments = FALSE)
  }
  if (Correction == "UnveilC" | Correction == "UnveiliC" | Correction == "UnveilJ") {
    entropy <- Shannon.ProbaVector(TunedPs, CheckArguments = FALSE)
    names(entropy) <- Correction
    return (entropy)
  }
  
  warning("Correction was not recognized")
  return (NA)
}
