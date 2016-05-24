# This function starts the estimation of the parameters. It runs for a fixed number of EM iterations.
# sigmaType:  Structure of the sigma matrices in the model.
# kKi:        Number of random effects per variance component.
# kLh:        Number of subvariance components in each variance component.
# kLhi:       Number of random effects in each subvariance component.
# kY, kX, kZ: Data and design matrices.
# EMit:       Number of EM iterations.
# MCit:       Number of intial MCMC iterations
# MCf:        Factor to increase the number of MCMC iterations.
# MCsd:       Standard deviation for the proposal step.

mcemMLE_t_fixed_df <- function (sigmaType, df, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM, controlTrust) {  
  # Number of fixed effects, random effects, variance and subvariance components.
  kP <- ncol(kX)
  kK <- ncol(kZ)
  kR <- length(kKi)
  u <- rep(0, kK)
  kL <- sum(kLh)
  
  # Parameters needed in sigma, one for diagonal, two for exchangeable and AR(1).
  beta <- initial[1:kP]
  sigma <- initial[-(1:kP)]
  
  QfunVal <- NULL
  theta <- c(beta, sigma)
  ovSigma <- constructSigma(pars = sigma, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
  outMLE <- matrix(0, controlEM$EMit, length(theta))
  outMLE[1, ] <- theta
  
  # MCMC step size tuning
  if (controlEM$MCsd == 0) {
    if (controlEM$verb >= 1)
      print("Tuning acceptance rate.")
    ar <- 1
    sdtune <- 1
    u <- rnorm(kK, rep(0, kK), sqrt(diag(ovSigma))) # Initial value for u
    while (ar > 0.4 | ar < 0.15) {
      uSample <- uSamplerCpp(beta = beta, sigma = ovSigma, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = 5000, sd0 = sdtune)
      ar <- length(unique(uSample[, 1])) / 5000
      if (ar < 0.15)
        sdtune <- 0.8 * sdtune
      if (ar > 0.4)
      sdtune <- 1.2 * sdtune
    }
    if (controlEM$verb >= 1)
      print(ar)
    controlEM$MCsd <- sdtune
  }

  # EM iterations
  j <- 2
  errorCounter <- 0
  while (j <= controlEM$EMit & sum(tail(errorCounter, 3)) < 3) {
    # Obtain MCMC sample for u with the current parameter estimates.
    u <- rnorm(kK, rep(0, kK), sqrt(diag(ovSigma))) # Initial value for u
    uSample <- uSamplerCpp(beta = beta, sigma = ovSigma, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = controlEM$MCit, sd0 = controlEM$MCsd)
    
    # Now we optimize.
    outTrust <- trust(toMaxDiag_t, parinit = theta, rinit = controlTrust$rinit, rmax = controlTrust$rmax, iterlim = controlTrust$iterlim, minimize = FALSE, u = uSample, sigmaType = sigmaType, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
    
    if (controlEM$verb >= 1)
      print(outTrust)
    outMLE[j, ] <- outTrust$argument
    QfunVal <- c(QfunVal, outTrust$value)
    
    # The current estimates are updated now
    beta <- outMLE[j, 1:kP]
    sigma <- outMLE[j, -c(1:kP)]
    theta <- c(beta, sigma)
    ovSigma <- constructSigma(pars = sigma, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
    if (controlEM$verb >= 1) {
      print(outMLE[1:j, ])
      if (controlEM$verb >= 2)
        print(ts.plot(uSample[, sample(1:kK, 1)]))
    }
    
    # Retuning the acceptance rate.
    ar <- length(unique(uSample[, 1]))/controlEM$MCit
    if (ar < 0.15 | ar > 0.4) {
      if (controlEM$verb >= 1)
        print("Tuning acceptance rate.")
      ar <- 1
      sdtune <- controlEM$MCsd
      u <- rnorm(kK, rep(0, kK), sqrt(diag(ovSigma))) # Initial value for u
      while (ar > 0.4 | ar < 0.15) {
        uSample.tmp <- uSamplerCpp(beta = beta, sigma = ovSigma, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = 5000, sd0 = sdtune)
        ar <- length(unique(uSample.tmp[, 1])) / 5000
        if (ar < 0.15)
          sdtune <- 0.8 * sdtune
        if (ar > 0.4)
          sdtune <- 1.2 * sdtune
      }
      if (controlEM$verb >= 1)
        print(ar)
      controlEM$MCsd <- sdtune
    }
    
    # Error checking
    error <- max(abs(outMLE[j, ] - outMLE[j - 1, ])/(abs(outMLE[j, ]) + controlEM$EMdelta))
    if(controlEM$verb >= 1)
      print(error)
    if (error < controlEM$EMepsilon) {
      errorCounter <- c(errorCounter, 1)
    } else {
      errorCounter <- c(errorCounter, 0)
    }
    
    # We modify the number of MCMC iterations
    if (j > 15 & controlEM$MCf < 1.1) {
      controlEM$MCf <- 1.2
    }
    if (sum(errorCounter) >= 2 | j > 30) {
      controlEM$MCf <- 1.5
    }
    controlEM$MCit <- controlEM$MCit * controlEM$MCf
    
    # Modify trust region
    controlTrust$rinit <- 2 * max(abs(outMLE[j, ] - outMLE[j - 1, ]))
    
    j <- j + 1
  }
  # Estimation of the information matrix.
  ovSigma <- constructSigma(pars = sigma, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
  B0 <- controlEM$MCit/controlEM$MCf
  uSample <- uSamplerCpp(beta = beta, sigma = ovSigma, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = B0, sd0 = controlEM$MCsd)
  iMatrix <- iMatrixDiagCpp_t(beta = beta, sigma = ovSigma, sigmaType = sigmaType, uSample = uSample, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = B0, sd0 = controlEM$MCsd)

  colnames(uSample) <- colnames(kZ)
  
  # loglikehood MCMC
  QfunMCMC <- MCMCloglikelihoodLogitCpp_t(beta = beta, sigma = ovSigma, sigmaType = sigmaType, u = uSample, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
  
  return(list(mcemEST = outMLE, iMatrix = iMatrix, QfunVal = QfunVal, QfunMCMC = QfunMCMC, randeff = uSample, y = kY, x = kX, z = kZ, EMerror = error, MCsd = controlEM$MCsd))
}