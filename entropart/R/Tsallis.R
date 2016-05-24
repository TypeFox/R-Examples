Tsallis <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("Tsallis")
}


Tsallis.ProbaVector <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
  dataTsallis <- -Ps^q * lnq(Ps, q)
  # Eliminate unobserved species
  dataTsallis[Ps == 0] <- 0
  
  entropy <- sum(dataTsallis)
  names(entropy) <- "None"
  return (entropy)
}


Tsallis.AbdVector <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return (bcTsallis(Ns=NorP, q=q, Correction=Correction, CheckArguments=CheckArguments))
}


Tsallis.integer <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcTsallis(Ns=NorP, q=q, Correction=Correction, CheckArguments=CheckArguments))
}


Tsallis.numeric <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
    return(Tsallis.ProbaVector(NorP, q=q, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return(Tsallis.AbdVector(NorP, q=q, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcTsallis <-
function(Ns, q = 1, Correction = "Best", CheckArguments = TRUE) 
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
    entropy <- Tsallis.ProbaVector(Ns/sum(Ns), q, CheckArguments = FALSE)
    names(entropy) <- Correction
    return (entropy)
  }
  
  
  # Common code for ZhangGrabchak. Useless if EntropyEstimation is used.
  # if (Correction == "ZhangGrabchak" | Correction == "ChaoWangJost") {
  #   Ps <- Ns/N
  #   V <- 1:(N-1)
  #   # p_V_Ns is an array, containing (1 - (n_s-1)/(n-j)) for each species (lines) and all j from 1 to n-1
  #   p_V_Ns <- outer(Ns, V, function(Ns, j) 1- (Ns-1)/(N-j))
  #   # Useful values are products from j=1 to v, so prepare cumulative products
  #   p_V_Ns <- t(apply(p_V_Ns, 1, cumprod))
  #   # Sum of products weighted by w_v
  #   S_v <- function(s) {
  #     Usedv <- 1:(N-Ns[s])
  #     return(sum(w_v[Usedv]*p_V_Ns[s, Usedv]))
  #   }
  # }
  
  # Shannon
  if (q == 1) {
    if (Correction == "Marcon") {
      ChaoShen <- bcShannon(Ns, Correction="ChaoShen", CheckArguments=FALSE)
      Grassberger <- bcShannon(Ns, Correction="Grassberger", CheckArguments=FALSE)
      entropy <- max(ChaoShen, Grassberger)
      names(entropy) <- Correction
      return (entropy)
    } else {
      if (Correction == "ZhangGrabchak") {
        # Weights. Useless if EntropyEstimation is used.
        # w_v <- 1/V
        # entropy <- sum(Ps*vapply(1:length(Ns), S_v, 0))
        entropy <- EntropyEstimation::Entropy.z(Ns)
        names(entropy) <- Correction
        return (entropy)
      } else {
      return (bcShannon(Ns, Correction=Correction, CheckArguments=FALSE)) 
      }
    }
  }
  
  # Not Shannon
  if (Correction == "ZhangGrabchak" | Correction == "ChaoWangJost") {
    # Weights. Useless here if EntropyEstimation is used, but weights are necessary for ChaoWangJost
    # i <- 1:N
    # w_vi <- (i-q)/i
    # w_v <- cumprod(w_vi)
    # ZhangGrabchak <- sum(Ps*vapply(1:length(Ns), S_v, 0))/(1-q)
    if (q==0) ZhangGrabchak <- length(Ns)-1 else ZhangGrabchak <- EntropyEstimation::Tsallis.z(Ns, q)
    # ZhangGrabchak stops here, but ChaoWangJost adds an estimator of the bias
    if (Correction == "ZhangGrabchak") {
      names(ZhangGrabchak) <- Correction
      return(ZhangGrabchak)
    } 
    # Calculate abundance distribution to obtain A
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
    # Eq 7d in Chao & Jost (2015). Terms for r in 1:(N-1) equal (-1)^r * w_v[r] * (A-1)^r. w_v is already available from ZhangGrabchak
    # Weights: here only if EntropyEstimation is used.
    i <- 1:N
    w_vi <- (i-q)/i
    w_v <- cumprod(w_vi)
    if (A == 1) {
      # The general formula of Eq 7d has a 0/0 part that must be forced to 0
      ChaoJostBias <- 0
    } else {
      Eq7dSum <- vapply(1:(N-1), function(r) w_v[r]*(1-A)^r, 0)
      # Calculate the estimator of the bias. Eq7dSum contains all terms of the sum except for r=0: the missing term equals 1.
      # The bias in Chao & Jost (2015) is that of the Hill number. It must be divided by 1-q to be applied to entropy.
      ChaoJostBias <- (Singletons/N*(1-A)^(1-N) * (A^(q-1) - sum(Eq7dSum) - 1))/(1-q)
    }
    entropy <- as.numeric(ZhangGrabchak + ChaoJostBias)
    names(entropy) <- Correction
    return (entropy)
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
    ChaoShen <- -sum(CPs^q * lnq(CPs, q) /(1 - (1-CPs)^N))
  }
  if (Correction == "ChaoShen" | Correction == "GenCov") {
    names(ChaoShen) <- Correction
    return(ChaoShen)
  }
  if (Correction == "Grassberger" | Correction == "Marcon") {
    Grassberger <- (1-N^(-q)*sum(Enq(Ns, q)))/(q-1)
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
  if (Correction == "Holste") {
    entropy <- 1/(1-q)*(beta(length(Ns)+N, q)*sum(1/beta(Ns+1, q))-1)
    names(entropy) <- Correction
    return (entropy)
  } 
  if (Correction == "Bonachela") {
    entropy <- 1/(1-q)*(beta(2+N, q)*sum(1/beta(Ns+1, q))-1)
    names(entropy) <- Correction
    return (entropy)
  }
  if (Correction == "UnveilC") {
    TunedPs <- as.ProbaVector(Ns, Correction="Chao2015", Unveiling="geom", RCorrection = "Chao1", CheckArguments=FALSE)
  }
  if (Correction == "UnveiliC") {
    TunedPs <- as.ProbaVector(Ns, Correction="Chao2015", Unveiling="geom", RCorrection = "iChao1", CheckArguments=FALSE)
  }
  if (Correction == "UnveilJ") {
    TunedPs <- as.ProbaVector(Ns, Correction="Chao2015", Unveiling="geom", RCorrection = "Jackknife", CheckArguments=FALSE)
  }
  if (Correction == "UnveilC" | Correction == "UnveiliC" | Correction == "UnveilJ") {
    entropy <- Tsallis.ProbaVector(TunedPs, q, CheckArguments = FALSE)
    names(entropy) <- Correction
    return (entropy)
  }
  
  warning("Correction was not recognized")
  return (NA)
}
