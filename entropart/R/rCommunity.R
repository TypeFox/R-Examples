rCommunity <- 
function(n, size = sum(NorP), NorP = 1, BootstrapMethod = "Chao2015",
         S = 300, Distribution = "lnorm", sd = 1, prob = 0.1, alpha=40, 
         CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  LSabundance <- function(N, alpha) {
    # Abundance of a species in a logseries distribution of size N and parameter alpha
    # adapted from Dan Lunn, http://www.stats.ox.ac.uk/~dlunn/BS1_05/BS1_Rcode.pdf
    x <- N/(N+alpha)
    u <- stats::runif(1)
    k <- 1
    P <- -(x)/log(1-x)
    F <- P
    while (u>=F) {
      P <- P*k*(x)/(k+1)
      k<-k+1
      F <- F+P
    }
    return(k)
  }
  
  
  Ps <- NULL
  # Draw probabilities
  if (length(NorP) == 1) {
    # Draw in a distribution
    Ps <- switch(Distribution,
                 geom = prob/(1-(1-prob)^S)*(1-prob)^(0:(S-1)),
                 lnorm = (stats::rlnorm(S, 0, sd) -> Ns)/sum(Ns),
                 lseries = (replicate(S, LSabundance(size, alpha))-> Ns)/sum(Ns),
                 bstick = c(cuts <- sort(stats::runif(S-1)), 1)- c(0, cuts)
    )
  } else {
    # Subsample
    if (abs(sum(NorP) - 1) < length(NorP)*.Machine$double.eps) {
      # Probabilities sum to 1, allowing rounding error
      Ps <- NorP    
    } else {
      Ns <- NorP
      # Abundances: Generate Ps according to the chosen method
      if (BootstrapMethod == "Chao2015") {
        Ps <- as.ProbaVector(Ns, Correction = "Chao2015", Unveiling = "geom", CheckArguments = FALSE)
      }
      if (BootstrapMethod == "Chao2013") {
        Ps <- as.ProbaVector(Ns, Correction = "Chao2013", Unveiling = "unif", CheckArguments = FALSE)
      }
      if (BootstrapMethod == "Marcon") {
        Ps <- Ns/sum(Ns)
      }
    }
  }
  
  if (is.null(Ps)) {
    warning ("The distribution to simulate has not been recognized")
    return(NA)
  }
  
  # Generate communities according to Ps
  Ns <- stats::rmultinom(n, size, Ps)
  if (n > 1) {
    # Return a MetaCommunity
    return(MetaCommunity(Ns))      
  } else {
    # Return a vector if a single community has been simulated
    return(as.AbdVector(Ns, Round = TRUE))
  }
}
