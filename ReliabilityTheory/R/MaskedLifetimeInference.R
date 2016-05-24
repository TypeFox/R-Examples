# t         = vector of system/network failure times
# signature = signature of the system/network for which inference is performed
#             May be a list of signatures which results in topology being inferred too
# cdfComp   = vectorised cumulative distribution function of component lifetime F_Y() with prototype:
#             function(y, parameters)
# pdfComp   = vectorised probability density function of component lifetime f_Y() with prototype:
#             function(y, parameters)
# rParmGivenData = generate random draws from f_{\Psi | Y} with prototype:
#                  function(y, ...)
# rCompGivenParm = generate random draws from f_{Y | \Psi} with prototype:
#                  function(parameters, t, censoring[-1,0,1])
# startParm = vector of starting values of named parameters in the correct order for the 'parameters' above
# iter      = number of MCMC iterations to perform
maskedInferenceIIDCustom <- function(t, signature, cdfComp, pdfComp, rParmGivenData, rCompGivenParm, startParm, iter, ...) {
  n <- length(t)
  # First off check if signature is a single sig/sig list/topology list (eg sccsO4)
  if( class(signature)=="list" && class(signature[[1]])=="list" && "signature"%in%names(signature[[1]]) && class(signature[[1]]$signature)=="numeric") {
    signature <- lapply(signature, `[[`, "signature")
  }
  sigIsList <- (class(signature)=="list")
  if(!sigIsList) { sig <- signature } # Here the signature selected randomly ("sig") never changes so assign off the bat
  numSigs <- length(signature)
  MgivenPsi <- rep(0, numSigs)
  
  # Setup return parameter matrix
  parm <- as.data.frame(matrix(startParm, nrow=iter+1, ncol=length(startParm), byrow=TRUE))
  names(parm) <- names(startParm)
  # and return topology vector
  M <- rep(NA, iter+1)
  
  progress <- 0
  for(it in 2:(iter+1)) {
    if(round(it*100/(iter+1), 0) != progress) { progress <- round(it*100/(iter+1), 0); cat("\r", progress, "% complete ...     "); }
    # DA block 1
    # DA block 1a
    # Sample topology given parameters and masked failure times
    if(sigIsList) {
      cdfT <- cdfComp(t, parm[it-1,,drop=FALSE], ...)
      pdfT <- pdfComp(t, parm[it-1,,drop=FALSE], ...)
      for(j in 1:numSigs) {
        i <- 1:length(signature[[j]])
        MgivenPsi[j] <- sum(vapply(1:n, function(ti) { log(sum(i * signature[[j]] * choose(length(signature[[j]]), i) * cdfT[ti]^(i-1) * (1-cdfT[ti])^(length(signature[[j]])-i) * pdfT[ti])) }, c(0)))
      }
      MgivenPsi <- MgivenPsi - max(MgivenPsi)
      M[it] <- sample(1:numSigs, 1, prob=exp(MgivenPsi))
      sig <- signature[[M[it]]]
    }
    
    # DA block 1b
    # Sample failure times of components given parameters
    y <- c()
    for(i in 1:n) {
      # Generate failure order statistics
      sig2 <- sig * dbinom(0:(length(sig)-1), length(sig)-1, cdfComp(t[i], parm[it-1,,drop=FALSE], ...))
      ##print(sig2)
      o <- sample.int(length(sig), size=1, prob=sig2)
      ##print(o)
      censoring <- c(rep(-1, o-1), 0, rep(1, length(sig)-o))
      y <- c(y, rCompGivenParm(parm[it-1,,drop=FALSE], t[i], censoring, ...))
    }
    ##cat("Simulated Components:", y, "\n")
    
    # DA block 2
    # Sample a new parameter values given simulated component failure times
    parm[it,] <- rParmGivenData(y, ...)
    ##cat("New Parameter Drawn:", parm[i], "\n===\n")
  }
  if(sigIsList) {
    return(cbind(topology=M[-1], parm[-1,,drop=FALSE]))
  }
  parm[-1,,drop=FALSE]
}


# t         = vector of n system/network failure times
# signature = signature of the system/network for which inference is performed
#             May be a list of signatures which results in topology being inferred too
# cdfComp   = vectorised cumulative distribution function of component lifetime F_Y() with prototype:
#             function(y, parameters)
# pdfComp   = vectorised probability density function of component lifetime f_Y() with prototype:
#             function(y, parameters)
# rParmGivenData = generate random draws from f_{\Psi | Y} with prototype:
#                  function(y, ...)
#                  Here the return must be a list: first element is like startHypPriorParm, second element is like startCompParm
# rCompGivenParm = generate random draws from f_{Y | \Psi} with prototype:
#                  function(parameters, t, censoring[-1,0,1])
# startCompParm = list of length n consisting of vectors of starting values of named parameters for component lifetime distribution in the correct order for the t and 'parameters' above
# startHypPriorParm = vector of starting values of named hyperprior parameters
# iter      = number of MCMC iterations to perform
maskedInferenceEXCHCustom <- function(t, signature, cdfComp, pdfComp, rParmGivenData, rCompGivenParm, startCompParm, startHypParm, iter, ...) {
  n <- length(t)
  # First off check if signature is a single sig/sig list/topology list (eg sccsO4)
  if( class(signature)=="list" && class(signature[[1]])=="list" && "signature"%in%names(signature[[1]]) && class(signature[[1]]$signature)=="numeric") {
    signature <- lapply(signature, `[[`, "signature")
  }
  sigIsList <- (class(signature)=="list")
  if(!sigIsList) { sig <- signature } # Here the signature selected randomly ("sig") never changes so assign off the bat
  numSigs <- length(signature)
  MgivenPsi <- rep(0, numSigs)
  
  # Setup parameter matrices
  psi <- lapply(1:n, function(i) { parm <- as.data.frame(matrix(startCompParm[[i]], nrow=iter+1, ncol=length(startCompParm[[i]]), byrow=TRUE)); names(parm) <- names(startCompParm[[i]]); parm })
  theta <- as.data.frame(matrix(startHypParm, nrow=iter+1, ncol=length(startHypParm), byrow=TRUE))
  names(theta) <- names(startHypParm)
  # and return topology vector
  M <- rep(NA, iter+1)
  
  progress <- 0
  for(it in 2:(iter+1)) {
    if(round(it*100/(iter+1), 0) != progress) { progress <- round(it*100/(iter+1), 0); cat("\r", progress, "% complete ...     "); }
    # DA block 1
    # DA block 1a
    # Sample topology given parameters and masked failure times
    if(sigIsList) {
      cdfT <- sapply(1:n, function(i) { cdfComp(t[i], psi[[i]][it-1,,drop=FALSE], ...) })
      pdfT <- sapply(1:n, function(i) { pdfComp(t[i], psi[[i]][it-1,,drop=FALSE], ...) })
      for(j in 1:numSigs) {
        i <- 1:length(signature[[j]])
        MgivenPsi[j] <- sum(vapply(1:n, function(ti) { log(sum(i * signature[[j]] * choose(length(signature[[j]]), i) * cdfT[ti]^(i-1) * (1-cdfT[ti])^(length(signature[[j]])-i) * pdfT[ti])) }, c(0)))
      }
      MgivenPsi <- MgivenPsi - max(MgivenPsi)
      M[it] <- sample(1:numSigs, 1, prob=exp(MgivenPsi))
      sig <- signature[[M[it]]]
    }
    
    # DA block 1b
    # Sample failure times of components given parameters
    m <- length(sig)
    y <- matrix(NA, nrow=n, ncol=m)
    for(i in 1:n) {
      # Generate failure order statistics
      sig2 <- exp(log(sig) + dbinom(0:(m-1), m-1, cdfComp(t[i], psi[[i]][it-1,,drop=FALSE], ...), log=TRUE))
      ##print(psi)
      o <- sample.int(m, size=1, prob=sig2)
      ##print(o)
      censoring <- c(rep(-1, o-1), 0, rep(1, m-o))
      y[i,] <- rCompGivenParm(psi[[i]][it-1,,drop=FALSE], t[i], censoring, ...)
    }
    ##cat("Simulated Components:", y, "\n")
    
    # DA block 2
    # Sample a new parameter values given simulated component failure times
    parm <- rParmGivenData(y, ...)
    theta[it,] <- parm[[1]]
    for(i in 1:n) {
      psi[[i]][it,] <- parm[[2]][[i]]
    }
    ##cat("New Parameter Drawn:", parm[i], "\n===\n")
  }
  if(sigIsList) {
    return(list(
      topology=M[-1],
      parameters=lapply(psi, function(psi_i) { psi_i[-1,,drop=FALSE]}),
      hyperparameters=theta[-1,,drop=FALSE]
    ))
  }
  return(list(
    parameters=lapply(psi, function(psi_i) { psi_i[-1,,drop=FALSE]}),
    hyperparameters=theta[-1,,drop=FALSE]
  ))
}
