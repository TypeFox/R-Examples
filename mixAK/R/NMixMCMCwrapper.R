##
##  PURPOSE:   (Reversible jump) MCMC for a normal mixture model
##             - wrapper to main simulation to allow vectorized call and parallel computation
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    06/11/2008
##              03/11/2011  parameter Cpar added to replace original z0, z1, censor, p, n, Cinteger, Cdouble
##              27/03/2015  mild revision to allow for factor covariates on mixture weights
##
##  FUNCTIONS:  NMixMCMCwrapper
##
## ======================================================================

## *************************************************************
## NMixMCMCwrapper
## *************************************************************
NMixMCMCwrapper <- function(chain = 1,
                            scale, prior, inits, Cpar, RJMCMC, CRJMCMC,
                            actionAll, nMCMC, keep.chains, PED,
                            dens.zero, lx_w)
{  
  thispackage <- "mixAK"

  p <- Cpar$dimy["p"]
  n <- Cpar$dimy["n"]
  nx_w <- Cpar$x_w[1]
  
  LTp <- p * (p + 1)/2
  Imat <- diag(p)
  rowsI <- row(Imat)[lower.tri(row(Imat), diag=TRUE)]
  colsI <- col(Imat)[lower.tri(col(Imat), diag=TRUE)] 
  naamLTp <- paste(".", rowsI, ".", colsI, sep="")      

  ##### scale:  derived variables
  ##### ============================================================================
  CshiftScale <- c(scale$shift, scale$scale)
  names(CshiftScale) <- c(paste("shift", 1:p, sep=""), paste("scale", 1:p, sep=""))

  ##### prior:  derived variables
  ##### ============================================================================
  CKmax <- as.numeric(prior$Kmax)

  ##### init:  derived variables
  ##### ============================================================================
  initz <- (inits[[chain]]$y - matrix(rep(scale$shift, n), ncol=p, byrow=TRUE))/matrix(rep(scale$scale, n), ncol=p, byrow=TRUE)
  
  CK <- as.numeric(inits[[chain]]$K)
  names(CK) <- "K"
  
  Cw <- c(inits[[chain]]$w, rep(0, CKmax - inits[[chain]]$K))
  if (nx_w == 1) names(Cw) <- paste("w", 1:CKmax, sep="") else names(Cw) <- paste(rep(paste("w", 1:CKmax, sep=""), nx_w), "-", rep(lx_w, each = CKmax), sep = "")
  
  if (p == 1){
    Cmu <- c(inits[[chain]]$mu, rep(0, CKmax - inits[[chain]]$K))
    names(Cmu) <- paste("mu", 1:CKmax, sep="")

    CLi <- c(inits[[chain]]$Li, rep(0, CKmax - inits[[chain]]$K))
    names(CLi) <- paste("Li", 1:CKmax, sep="")        
    
  }else{
    Cmu <- c(t(inits[[chain]]$mu), rep(0, p*(CKmax - inits[[chain]]$K)))
    names(Cmu) <- paste("mu", rep(1:CKmax, each=p), ".", rep(1:p, CKmax), sep="")

    CLi <- c(inits[[chain]]$Li, rep(0, LTp*(CKmax - inits[[chain]]$K)))
    names(CLi) <- paste("Li", rep(1:CKmax, each=LTp), rep(naamLTp, CKmax), sep="")
  }  

  CgammaInv <- inits[[chain]]$gammaInv

  Cr <- inits[[chain]]$r - 1
 
  ##### Some additional parameters
  ##### =============================================================================
  if (prior$priorK == "fixed") lsum_Ir <- n * CK
  else                         lsum_Ir <- 1
  
  ########## ========== MCMC sampling ========== ##########
  ########## =================================== ##########
  cat(paste("\nChain number ", chain, "\n==============\n", sep=""))
  cat(paste("MCMC sampling started on ", date(), ".\n", sep=""))
  MCMC <- .C("NMix_MCMC",
             z0                   = as.double(t(Cpar$z0)),
             z1                   = as.double(t(Cpar$z1)),
             censor               = as.integer(t(Cpar$censor)),
             nxw_xw               = as.integer(Cpar$x_w),
             dimy                 = as.integer(Cpar$dimy),
             shiftScale           = as.double(CshiftScale),             
             nMCMC                = as.integer(nMCMC),
             priorInt             = as.integer(Cpar$priorInt),
             priorDouble          = as.double(Cpar$priorDouble),
             priorRJMCMC          = as.double(CRJMCMC),
             priorRJMCMCint       = as.integer(actionAll),
             z                    = as.double(t(initz)),
             z_first              = double(p*n),
             K                    = as.integer(CK),             
             w                    = as.double(Cw),
             mu                   = as.double(Cmu),
             Q                    = double(CKmax * LTp),
             Sigma                = double(CKmax * LTp),
             Li                   = as.double(CLi),
             gammaInv             = as.double(CgammaInv),
             r                    = as.integer(Cr),
             r_first              = integer(n),
             chK                  = integer(nMCMC["keep"]),
             chw                  = double(CKmax * nMCMC["keep"] * nx_w),
             chmu                 = double(p * CKmax * nMCMC["keep"]),
             chQ                  = double(LTp * CKmax * nMCMC["keep"]),
             chSigma              = double(LTp * CKmax * nMCMC["keep"]),             
             chLi                 = double(LTp * CKmax * nMCMC["keep"]),
             chgammaInv           = double(p * nMCMC["keep"]),
             chorder              = integer(CKmax * nMCMC["keep"]),
             chrank               = integer(CKmax * nMCMC["keep"]),             
             chMean               = double(p * nMCMC["keep"] * nx_w),
             chCorr               = double(LTp * nMCMC["keep"] * nx_w),
             chMeanData           = double(p * nMCMC["keep"] * nx_w),             
             chCorrData           = double(LTp * nMCMC["keep"] * nx_w),
             chLogL0              = double(nMCMC["keep"]),
             chLogL1              = double(nMCMC["keep"]),
             chDevCompl           = double(nMCMC["keep"]),             
             chDevObs             = double(nMCMC["keep"]),
             chDevCompl.inHat     = double(nMCMC["keep"]),
             pm.z                 = double(p * n),
             pm.indLogL0          = double(n),
             pm.indLogL1          = double(n),
             pm.indDevCompl       = double(n),             
             pm.indDevObs         = double(n),
             pm.indDevCompl.inHat = double(n),             
             pm.pred.dens         = double(n),             
             pm.w                 = double(CKmax * nx_w),
             pm.mu                = double(p * CKmax),
             pm.Q                 = double(LTp * CKmax),
             pm.Sigma             = double(LTp * CKmax),
             pm.Li                = double(LTp * CKmax),
             sum_Ir               = integer(lsum_Ir),
             sum_Pr_y             = double(lsum_Ir),
             iter                 = as.integer(0),
             nMoveAccept          = integer(9),
             err                  = as.integer(0),
             PACKAGE=thispackage)
  cat(paste("MCMC sampling finished on ", date(), ".\n", sep=""))
  if (MCMC$err) stop("Something went wrong.")

  
  ########## ========== State of MCMC (last and first kept) ========== ##########
  ########## ========================================================= ##########  
  if (p == 1){
    state.mu       <- as.numeric(MCMC$mu[1:MCMC$K])
    state_first.mu <- as.numeric(MCMC$chmu[1:MCMC$chK[1]])    
    names(state.mu)       <- paste("mu", 1:MCMC$K, sep="")
    names(state_first.mu) <- paste("mu", 1:MCMC$chK[1], sep="")    
    
    state.Li       <- as.numeric(MCMC$Li[1:MCMC$K])
    state_first.Li <- as.numeric(MCMC$chLi[1:MCMC$chK[1]])    
    names(state.Li)       <- paste("Li", 1:MCMC$K, sep="")
    names(state_first.Li) <- paste("Li", 1:MCMC$chK[1], sep="")
    
    state.Sigma       <- (1 / state.Li)^2
    state_first.Sigma <- (1 / state_first.Li)^2
    names(state.Sigma)       <- paste("Sigma", 1:MCMC$K, sep="")
    names(state_first.Sigma) <- paste("Sigma", 1:MCMC$chK[1], sep="")
    
    state.Q       <- as.numeric(MCMC$Q[1:MCMC$K])
    state_first.Q <- as.numeric(MCMC$chQ[1:MCMC$chK[1]])    
    names(state.Q)       <- paste("Q", 1:MCMC$K, sep="")
    names(state_first.Q) <- paste("Q", 1:MCMC$chK[1], sep="")        
  }else{
    state.mu       <- matrix(MCMC$mu[1:(p*MCMC$K)], ncol=p, byrow=TRUE)
    state_first.mu <- matrix(MCMC$chmu[1:(p*MCMC$chK[1])], ncol=p, byrow=TRUE)    
    rownames(state.mu)       <- paste("j", 1:MCMC$K, sep="")
    rownames(state_first.mu) <- paste("j", 1:MCMC$chK[1], sep="")    
    colnames(state.mu) <- colnames(state_first.mu) <- paste("m", 1:p, sep="")
    
    state.Li       <- as.numeric(MCMC$Li[1:(LTp*MCMC$K)])
    state_first.Li <- as.numeric(MCMC$chLi[1:(LTp*MCMC$chK[1])])    
    names(state.Li)       <- paste("Li", rep(1:MCMC$K, each=LTp), rep(naamLTp, MCMC$K), sep="")
    names(state_first.Li) <- paste("Li", rep(1:MCMC$chK[1], each=LTp), rep(naamLTp, MCMC$chK[1]), sep="")    
    
    state.Sigma <- matrix(NA, ncol=p, nrow=p*MCMC$K)
    rownames(state.Sigma) <- paste("j", rep(1:MCMC$K, each=p), ".", rep(1:p, MCMC$K), sep="")
    colnames(state.Sigma) <- paste("m", 1:p, sep="")            
    for (j in 1:MCMC$K){
      tmpSigma <- matrix(0, nrow=p, ncol=p)
      tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- state.Li[((j-1)*LTp+1):(j*LTp)]
      tmpSigma <- tmpSigma %*% t(tmpSigma)
      tmpSigma <- chol2inv(chol(tmpSigma))
      state.Sigma[((j-1)*p+1):(j*p),] <- tmpSigma
    }

    state_first.Sigma <- matrix(NA, ncol=p, nrow=p*MCMC$chK[1])
    rownames(state_first.Sigma) <- paste("j", rep(1:MCMC$chK[1], each=p), ".", rep(1:p, MCMC$chK[1]), sep="")
    colnames(state_first.Sigma) <- paste("m", 1:p, sep="")            
    for (j in 1:MCMC$chK[1]){
      tmpSigma <- matrix(0, nrow=p, ncol=p)
      tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- state_first.Li[((j-1)*LTp+1):(j*LTp)]
      tmpSigma <- tmpSigma %*% t(tmpSigma)
      tmpSigma <- chol2inv(chol(tmpSigma))
      state_first.Sigma[((j-1)*p+1):(j*p),] <- tmpSigma
    }
    
    state.Q       <- as.numeric(MCMC$Q[1:(LTp*MCMC$K)])
    state_first.Q <- as.numeric(MCMC$chQ[1:(LTp*MCMC$chK[1])])    
    names(state.Q)       <- paste("Q", rep(1:MCMC$K, each=LTp), rep(naamLTp, MCMC$K), sep="")
    names(state_first.Q) <- paste("Q", rep(1:MCMC$chK[1], each=LTp), rep(naamLTp, MCMC$chK[1]), sep="")        
  }
  nCompTotal <- sum(MCMC$chK)
  freqK <- table(MCMC$chK)
  propK <- prop.table(freqK)

  state.z       <- matrix(MCMC$z, ncol=p, byrow=TRUE)
  state_first.z <- matrix(MCMC$z_first, ncol=p, byrow=TRUE)  
  state.y       <- state.z       * matrix(rep(scale$scale, n), ncol=p, byrow=TRUE) + matrix(rep(scale$shift, n), ncol=p, byrow=TRUE)
  state_first.y <- state_first.z * matrix(rep(scale$scale, n), ncol=p, byrow=TRUE) + matrix(rep(scale$shift, n), ncol=p, byrow=TRUE)  
  
  ########## ========== Deviance and DIC's ========== ##########
  ########## ======================================== ##########
     ### when computing Dbar, truncate Inf values
  detS      <- prod(scale$scale)
  idetS     <- 1 / detS
  log.idetS <- -log(detS)

  REMOVE <- MCMC$chDevObs >= Inf
  D3bar   <- mean(MCMC$chDevObs[!REMOVE]) - 2*n*log.idetS
  D3inBar <- -2 * sum(log(MCMC$pm.pred.dens)) - 2*n*log.idetS
  pD3     <- D3bar - D3inBar
  DIC3    <- D3bar + pD3

## 06/03/2008:  Do not return DIC4 as it is not clear to me whether it is correct  
#  if (prior$priormuQ == "naturalC"){  
#    D4bar   <- mean(MCMC$chDevCompl[MCMC$chDevCompl < Inf]) - 2*n*log.idetS
#    D4inBar <- mean(MCMC$chDevCompl.inHat) - 2*n*log.idetS
#    pD4     <- D4bar - D4inBar
#    DIC4    <- D4bar + pD4
#    DIC <- data.frame(DIC=c(DIC3, DIC4), pD=c(pD3, pD4), D.bar=c(D3bar, D4bar), D.in.bar=c(D3inBar, D4inBar))
#    rownames(DIC) <- paste("DIC", 3:4, sep="")
#  }else{
    DIC <- data.frame(DIC=DIC3, pD=pD3, D.bar=D3bar, D.in.bar=D3inBar)
    rownames(DIC) <- paste("DIC", 3, sep="")    
#  }  

    
  ########## ========== nMoveAccept ============ ##########
  ########## =================================== ##########  
  names(MCMC$nMoveAccept) <- c("nGibbsK", "nSplit", "nCombine", "nBirth", "nDeath", "nAcceptSplit", "nAcceptCombine", "nAcceptBirth", "nAcceptDeath")
  if (prior$priorK == "fixed"){
    nMove <- MCMC$nMoveAccept["nGibbsK"]
    nAccept <- MCMC$nMoveAccept["nGibbsK"]
    percAccept <- 100
    moves <- data.frame(nMove, nAccept, percAccept)
    rownames(moves) <- "Gibbs with fixed K"
  }else{
    nMove <- MCMC$nMoveAccept[1:5]
    nAccept <- c(MCMC$nMoveAccept[1], MCMC$nMoveAccept[6:9])
    percAccept <- round((nAccept/nMove)*100, 2)
    moves <- data.frame(nMove, nAccept, percAccept)
    rownames(moves) <- c("Gibbs with fixed K", "Split", "Combine", "Birth", "Death")
  }  
  colnames(moves) <- c("Performed", "Accepted", "Proportion accepted (%)")
    
  ########## ========== Put everything together ========== ##########
  ########## ============================================= ##########      
  RET <- list(iter          = MCMC$iter,
              nMCMC         = nMCMC,
              dim           = p,
              nx_w          = nx_w,
              lx_w          = lx_w,
              #x_w           = if (nx_w > 1) Cpar$x_w[-1] else 0,
              prior         = prior,
              init          = inits[[chain]],
              state.first   = list(y        = state_first.y,
                                   K        = as.numeric(MCMC$chK[1]),
                                   w        = as.numeric(MCMC$chw[1:(MCMC$chK[1] * nx_w)]),
                                   mu       = state_first.mu,
                                   Li       = state_first.Li,
                                   Sigma    = state_first.Sigma,
                                   Q        = state_first.Q,                
                                   gammaInv = as.numeric(MCMC$chgammaInv[1:p]),
                                   r        = as.numeric(MCMC$r_first+ 1)),              
              state.last    = list(y        = state.y,
                                   K        = as.numeric(MCMC$K),
                                   w        = as.numeric(MCMC$w[1:(MCMC$K * nx_w)]),
                                   mu       = state.mu,
                                   Li       = state.Li,
                                   Sigma    = state.Sigma,
                                   Q        = state.Q,                
                                   gammaInv = as.numeric(MCMC$gammaInv),
                                   r        = as.numeric(MCMC$r + 1)),              
              RJMCMC        = RJMCMC,
              scale         = scale,
              freqK         = freqK,
              propK         = propK,
              DIC           = DIC,
              moves         = moves)
  if (nx_w > 1){
    names(RET$state.last$w)  <- paste("w", rep(1:MCMC$K, nx_w), "-", rep(lx_w, each = MCMC$K), sep="")
    names(RET$state.first$w) <- paste("w", rep(1:MCMC$chK[1], nx_w), "-", rep(lx_w, each = MCMC$chK[1]), sep="")
  }else{
    names(RET$state.last$w)  <- paste("w", 1:MCMC$K, sep="")
    names(RET$state.first$w) <- paste("w", 1:MCMC$chK[1], sep="")  
  }    
  names(RET$state.last$gammaInv) <- names(RET$state.first$gammaInv) <- paste("gammaInv", 1:p, sep="")
  names(RET$state.last$r)        <- names(RET$state.first$r)        <- paste("r", 1:n, sep="")

  
  ########## ========== Chains for basic parameters ========== ##########
  ########## ================================================= ##########
  if (keep.chains | PED){
    RET$K <- as.numeric(MCMC$chK)
    MCMC$chK <- NULL
  
    RET$w <- as.numeric(MCMC$chw[1:(nx_w * nCompTotal)])
    MCMC$chw <- NULL

    RET$mu <- as.numeric(MCMC$chmu[1:(p*nCompTotal)])
    MCMC$chmu <- NULL

    RET$Li <- as.numeric(MCMC$chLi[1:(LTp*nCompTotal)])
    MCMC$chLi <- NULL    
  }

  ########## ========== Chains for mixture (overall) means, std. deviations and correlations ========== ##########
  ########## ========================================================================================== ##########    
  MCMC$chMean     <- matrix(MCMC$chMean,     ncol = p   * nx_w, byrow = TRUE)
  MCMC$chCorr     <- matrix(MCMC$chCorr,     ncol = LTp * nx_w, byrow = TRUE)
  MCMC$chMeanData <- matrix(MCMC$chMeanData, ncol = p   * nx_w, byrow = TRUE)
  MCMC$chCorrData <- matrix(MCMC$chCorrData, ncol = LTp * nx_w, byrow = TRUE)
  if (nx_w > 1){
    colnames(MCMC$chMean) <- paste("z.Mean.", rep(1:p, nx_w), "-", rep(lx_w, each = p), sep="")
    colnames(MCMC$chCorr) <- paste("z.Corr", rep(naamLTp, nx_w), "-", rep(lx_w, each = LTp), sep="")
    colnames(MCMC$chMeanData) <- paste("y.Mean.", rep(1:p, nx_w), "-", rep(lx_w, each = p), sep="")
    colnames(MCMC$chCorrData) <- paste("y.Corr", rep(naamLTp, nx_w), "-", rep(lx_w, each = LTp), sep="")
    for (ixw in 0:(nx_w - 1)){
      colnames(MCMC$chCorr)[ixw * LTp + ((0:(p-1))*(2*p - (0:(p-1)) + 1))/2 + 1]     <- paste("z.SD.", 1:p, "-", lx_w[ixw + 1], sep="")
      colnames(MCMC$chCorrData)[ixw * LTp + ((0:(p-1))*(2*p - (0:(p-1)) + 1))/2 + 1] <- paste("y.SD.", 1:p, "-", lx_w[ixw + 1], sep="")      
    }    
  }else{    
    colnames(MCMC$chMean) <- paste("z.Mean.", 1:p, sep="")
    colnames(MCMC$chCorr) <- paste("z.Corr", naamLTp, sep="")
    colnames(MCMC$chCorr)[((0:(p-1))*(2*p - (0:(p-1)) + 1))/2 + 1] <- paste("z.SD.", 1:p, sep="")
    colnames(MCMC$chMeanData) <- paste("y.Mean.", 1:p, sep="")
    colnames(MCMC$chCorrData) <- paste("y.Corr", naamLTp, sep="")
    colnames(MCMC$chCorrData)[((0:(p-1))*(2*p - (0:(p-1)) + 1))/2 + 1] <- paste("y.SD.", 1:p, sep="")
  }  

  
  ######### =========== Work-out the chains further ========== ##########
  ######### ================================================== ##########
  if (keep.chains){
    RET$mixture <- as.data.frame(cbind(MCMC$chMeanData, MCMC$chCorrData, MCMC$chMean, MCMC$chCorr))
    
    RET$Q <- as.numeric(MCMC$chQ[1:(LTp*nCompTotal)])
    MCMC$chQ <- NULL

    RET$Sigma <- as.numeric(MCMC$chSigma[1:(LTp*nCompTotal)])
    MCMC$chSigma <- NULL

    RET$gammaInv <- matrix(MCMC$chgammaInv, ncol=p, byrow=TRUE)
    colnames(RET$gammaInv) <- paste("gammaInv", 1:p, sep="")  
    MCMC$chgammaInv <- NULL
  
    RET$order <- as.numeric(MCMC$chorder[1:nCompTotal] + 1)
    MCMC$chorder <- NULL

    RET$rank <- as.numeric(MCMC$chrank[1:nCompTotal] + 1)
    MCMC$chrank <- NULL
  
    ########## ========== Chains for deviances and related quantities  ========== ##########
    ########## ================================================================== ##########
  ## 06/03/2008:  Do not return DIC4 related quantities as it is not clear to me whether it is correct    
  #  if (prior$priormuQ == "naturalC"){
  #    RET$deviance <- data.frame(LogL0              = MCMC$chLogL0 + log.idetS,
  #                               LogL1              = MCMC$chLogL1,
  #                               dev.complete       = MCMC$chDevCompl - 2*n*log.idetS,
  #                               dev.observed       = MCMC$chDevObs - 2*n*log.idetS,
  #                               dev.complete.inHat = MCMC$chDevCompl.inHat - 2*n*log.idetS)
  #  }else{
      RET$deviance <- data.frame(LogL0        = MCMC$chLogL0 + log.idetS,
                                 LogL1        = MCMC$chLogL1,
                                 dev.complete = MCMC$chDevCompl - 2*n*log.idetS,
                                 dev.observed = MCMC$chDevObs - 2*n*log.idetS)
  #  }
    MCMC$chLogL0          <- NULL
    MCMC$chLogL1          <- NULL
    MCMC$chDevCompl       <- NULL
    MCMC$chDevObs         <- NULL
    MCMC$chDevCompl.inHat <- NULL
  }else{
    MCMC$chK <- NULL
    MCMC$chw <- NULL
    MCMC$chmu <- NULL
    MCMC$chQ <- NULL
    MCMC$chSigma <- NULL
    MCMC$chLi <- NULL
    MCMC$chgammaInv <- NULL
    MCMC$chorder <- NULL
    MCMC$chrank <- NULL

    MCMC$chLogL0          <- NULL
    MCMC$chLogL1          <- NULL
    MCMC$chDevCompl       <- NULL
    MCMC$chDevObs         <- NULL
    MCMC$chDevCompl.inHat <- NULL        
  }  

  ########## ========== Clustering based on posterior P(alloc = k | y) or on P(alloc = k | theta, y)       ========== ##########
  ########## ======================================================================================================== ##########
  if (prior$priorK == "fixed"){
    if (CK == 1){
      RET$poster.comp.prob_u <- RET$poster.comp.prob_b <- matrix(1, nrow = n, ncol = 1)
    }else{

      ### Using mean(I(r=k))
      MCMC$sum_Ir <- matrix(MCMC$sum_Ir, ncol = CK, nrow = n, byrow = TRUE)
      Denom <- apply(MCMC$sum_Ir, 1, sum)
      RET$poster.comp.prob_u <- MCMC$sum_Ir / matrix(rep(Denom, CK), ncol = CK, nrow = n)

      ### Using mean(P(r=k | theta, b, y))
      MCMC$sum_Pr_y <- matrix(MCMC$sum_Pr_y, ncol = CK, nrow = n, byrow = TRUE)
      RET$poster.comp.prob_b <- MCMC$sum_Pr_y/ matrix(rep(Denom, CK), ncol = CK, nrow = n)        
    }  
  }  
  
  ########## ========== Posterior means of scaled and original observations (useful in the case of censoring)  ========== ##########
  ########## ============================================================================================================ ##########  
  MCMC$pm.z <- matrix(MCMC$pm.z, ncol=p, byrow=TRUE)
  RET$pm.y <- MCMC$pm.z * matrix(rep(scale$scale, n), ncol=p, byrow=TRUE) + matrix(rep(scale$shift, n), ncol=p, byrow=TRUE)
  RET$pm.y <- as.data.frame(RET$pm.y)    
  RET$pm.z <- as.data.frame(MCMC$pm.z)
  colnames(RET$pm.y) <- paste("y", 1:p, sep="")    
  colnames(RET$pm.z) <- paste("z", 1:p, sep="")
  MCMC$pm.z <- NULL
  
  ########## ========== Posterior means of deviance contributions  ========== ##########
  ########## ================================================================ ##########
## 06/03/2008:  Do not return DIC4 related quantities as it is not clear to me whether it is correct      
#  if (prior$priormuQ == "naturalC"){
#    RET$pm.indDev <- data.frame(LogL0              = MCMC$pm.indLogL0 + log.idetS,
#                                LogL1              = MCMC$pm.indLogL1,
#                                dev.complete       = MCMC$pm.indDevCompl - 2*log.idetS,
#                                dev.observed       = MCMC$pm.indDevObs - 2*log.idetS,
#                                dev.complete.inHat = MCMC$pm.indDevCompl.inHat - 2*log.idetS)
#  }else{
    RET$pm.indDev <- data.frame(LogL0        = MCMC$pm.indLogL0 + log.idetS,
                                LogL1        = MCMC$pm.indLogL1,
                                dev.complete = MCMC$pm.indDevCompl - 2*log.idetS,
                                dev.observed = MCMC$pm.indDevObs - 2*log.idetS)
#  }
  MCMC$pm.LogL0             <- NULL
  MCMC$pm.LogL1             <- NULL
  MCMC$pm.indDevCompl       <- NULL
  MCMC$pm.indDevObs         <- NULL
  MCMC$pm.indDevCompl.inHat <- NULL

  ########## ========== Predictive density evaluated at observed data =============== ##########
  ########## ======================================================================== ##########
  RET$pred.dens <- idetS * MCMC$pm.pred.dens
  MCMC$pm.pred.dens <- NULL  
  
  ########## ========== Some adjustment of the output when priorK is fixed ========== ##########
  ########## ======================================================================== ##########    
  if (prior$priorK == "fixed"){                              ### Chains for order, rank, w, mu, Q, Sigma, Li are returned in a form of matrices
    if (keep.chains | PED){
      RET$w <- matrix(RET$w, ncol = CKmax * nx_w, byrow = TRUE)
      if (nx_w > 1){
        colnames(RET$w) <- paste("w", rep(1:CKmax, nx_w), "-", rep(lx_w, each = CKmax), sep = "")
      }else{
        colnames(RET$w) <- paste("w", 1:CKmax, sep = "")
      }    

      RET$mu <- matrix(RET$mu, ncol = p * CKmax, byrow = TRUE)
      colnames(RET$mu) <- paste("mu.", rep(1:CKmax, each = p), ".", rep(1:p, CKmax), sep = "")
    
      RET$Li <- matrix(RET$Li, ncol = LTp*CKmax, byrow = TRUE)
      colnames(RET$Li) <- paste("Li", rep(1:CKmax, each = LTp), rep(naamLTp, CKmax), sep = "")
    }  
    
    if (keep.chains){
      RET$order <- matrix(RET$order, ncol=CKmax, byrow=TRUE)
      colnames(RET$order) <- paste("order", 1:CKmax, sep="")

      RET$rank <- matrix(RET$rank, ncol=CKmax, byrow=TRUE)
      colnames(RET$rank) <- paste("rank", 1:CKmax, sep="")

      RET$Q <- matrix(RET$Q, ncol=LTp*CKmax, byrow=TRUE)
      colnames(RET$Q) <- paste("Q", rep(1:CKmax, each=LTp), rep(naamLTp, CKmax), sep="")

      RET$Sigma <- matrix(RET$Sigma, ncol=LTp*CKmax, byrow=TRUE)
      colnames(RET$Sigma) <- paste("Sigma", rep(1:CKmax, each=LTp), rep(naamLTp, CKmax), sep="")      
    }
  }  
  
  ########## ========== Posterior summary statistics for moments of the scaled and unscaled mixture ========== ##########
  ########## ================================================================================================= ##########
  qProbs <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
  nSumm <- c("Mean", "Std.Dev.", "Min.", "2.5%", "1st Qu.", "Median", "3rd Qu.", "97.5%", "Max.")
  if (p == 1 & nx_w == 1){
    meanmix.Mean  <- mean(MCMC$chMean, na.rm=TRUE)
    quantmix.Mean <- quantile(MCMC$chMean, prob=qProbs, na.rm=TRUE)    
    sdmix.Mean   <- sd(MCMC$chMean, na.rm=TRUE)
    smix.Mean <- c(meanmix.Mean, sdmix.Mean, quantmix.Mean)
    names(smix.Mean) <- nSumm

    meanmix.SDCorr  <- mean(MCMC$chCorr, na.rm=TRUE)
    quantmix.SDCorr <- quantile(MCMC$chCorr, prob=qProbs, na.rm=TRUE)    
    sdmix.SDCorr   <- sd(MCMC$chCorr, na.rm=TRUE)
    smix.SDCorr <- c(meanmix.SDCorr, sdmix.SDCorr, quantmix.SDCorr)
    names(smix.SDCorr) <- nSumm

    meany.Mean  <- mean(MCMC$chMeanData, na.rm=TRUE)
    quanty.Mean <- quantile(MCMC$chMeanData, prob=qProbs, na.rm=TRUE)    
    sdy.Mean   <- sd(MCMC$chMeanData, na.rm=TRUE)
    sy.Mean <- c(meany.Mean, sdy.Mean, quanty.Mean)
    names(sy.Mean) <- nSumm

    meany.SDCorr  <- mean(MCMC$chCorrData, na.rm=TRUE)
    quanty.SDCorr <- quantile(MCMC$chCorrData, prob=qProbs, na.rm=TRUE)    
    sdy.SDCorr   <- sd(MCMC$chCorrData, na.rm=TRUE)
    sy.SDCorr <- c(meany.SDCorr, sdy.SDCorr, quanty.SDCorr)
    names(sy.SDCorr) <- nSumm    
  }else{
    meanmix.Mean <- apply(MCMC$chMean, 2, mean, na.rm=TRUE)
    quantmix.Mean <- apply(MCMC$chMean, 2, quantile, prob=qProbs, na.rm=TRUE)    
    sdmix.Mean <- apply(MCMC$chMean, 2, sd, na.rm=TRUE)
    smix.Mean <- rbind(meanmix.Mean, sdmix.Mean, quantmix.Mean)
    rownames(smix.Mean) <- nSumm

    meanmix.SDCorr <- apply(MCMC$chCorr, 2, mean, na.rm=TRUE)
    quantmix.SDCorr <- apply(MCMC$chCorr, 2, quantile, prob=qProbs, na.rm=TRUE)    
    sdmix.SDCorr <- apply(MCMC$chCorr, 2, sd, na.rm=TRUE)
    smix.SDCorr <- rbind(meanmix.SDCorr, sdmix.SDCorr, quantmix.SDCorr)
    rownames(smix.SDCorr) <- nSumm

    meany.Mean <- apply(MCMC$chMeanData, 2, mean, na.rm=TRUE)
    quanty.Mean <- apply(MCMC$chMeanData, 2, quantile, prob=qProbs, na.rm=TRUE)    
    sdy.Mean <- apply(MCMC$chMeanData, 2, sd, na.rm=TRUE)
    sy.Mean <- rbind(meany.Mean, sdy.Mean, quanty.Mean)
    rownames(sy.Mean) <- nSumm

    meany.SDCorr <- apply(MCMC$chCorrData, 2, mean, na.rm=TRUE)
    quanty.SDCorr <- apply(MCMC$chCorrData, 2, quantile, prob=qProbs, na.rm=TRUE)    
    sdy.SDCorr <- apply(MCMC$chCorrData, 2, sd, na.rm=TRUE)
    sy.SDCorr <- rbind(meany.SDCorr, sdy.SDCorr, quanty.SDCorr)
    rownames(sy.SDCorr) <- nSumm        
  }  
  
  RET$summ.y.Mean <- sy.Mean
  RET$summ.y.SDCorr <- sy.SDCorr  
  
  RET$summ.z.Mean <- smix.Mean
  RET$summ.z.SDCorr <- smix.SDCorr

  MCMC$chMeanData <- NULL
  MCMC$chCorrData <- NULL
  MCMC$chMean <- NULL
  MCMC$chCorr <- NULL 
  
  ########## ========== Posterior means for mixture components when priorK == "fixed" ========== ##########
  ########## =================================================================================== ##########  
  if (prior$priorK == "fixed"){
                                                      ##### I am not sure whether the posterior means (especially of variance components) are useful!
                                                      ##### In any case, they should be used with care
                                                      ##### -----------------------------------------------------------------------------------------
    RET$poster.mean.w <- as.numeric(MCMC$pm.w)
    if (nx_w > 1){
      names(RET$poster.mean.w) <- paste("w", rep(1:CKmax, nx_w), "-", rep(lx_w, each = CKmax), sep = "")
    }else{
      names(RET$poster.mean.w) <- paste("w", 1:CKmax, sep="")
    }    

    RET$poster.mean.mu <- matrix(MCMC$pm.mu, nrow=CKmax, ncol=p, byrow=TRUE)
    rownames(RET$poster.mean.mu) <- paste("j", 1:CKmax, sep="")
    colnames(RET$poster.mean.mu) <- paste("m", 1:p, sep="")

    RET$poster.mean.Q <- RET$poster.mean.Sigma <- RET$poster.mean.Li <- list()
    for (j in 1:CKmax){
      tmpQ <- matrix(0, nrow=p, ncol=p)
      tmpQ[lower.tri(tmpQ, diag=TRUE)] <- MCMC$pm.Q[((j-1)*LTp+1):(j*LTp)]
      tmpQ[upper.tri(tmpQ, diag=FALSE)] <- t(tmpQ)[upper.tri(t(tmpQ), diag=FALSE)]
      RET$poster.mean.Q[[j]] <- tmpQ
      
      tmpSigma <- matrix(0, nrow=p, ncol=p)
      tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- MCMC$pm.Sigma[((j-1)*LTp+1):(j*LTp)]
      tmpSigma[upper.tri(tmpSigma, diag=FALSE)] <- t(tmpSigma)[upper.tri(t(tmpSigma), diag=FALSE)]
      RET$poster.mean.Sigma[[j]] <- tmpSigma
      
      tmpLi <- matrix(0, nrow=p, ncol=p)
      tmpLi[lower.tri(tmpLi, diag=TRUE)] <- MCMC$pm.Li[((j-1)*LTp+1):(j*LTp)]
      RET$poster.mean.Li[[j]] <- tmpLi      
    }
    names(RET$poster.mean.Q) <- names(RET$poster.mean.Sigma) <- names(RET$poster.mean.Li) <- paste("j", 1:CKmax, sep="")    
  }  

  ########## ========== Additional objects (added on 08/02/2010) ========== ##########
  ########## ============================================================== ##########
  RET$relabel <- list(type = "mean", par = 1)       #### default re-labeling is performed using the first margin of the mixture means
  RET$Cpar <- Cpar
  
  class(RET) <- "NMixMCMC"
  return(RET)
}  

