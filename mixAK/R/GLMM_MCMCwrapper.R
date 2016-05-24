##
##  PURPOSE:   Generalized linear mixed model with possibly several response variables
##             and normal mixtures in the distribution of random effects
##             - wrapper to main simulation to allow vectorized call and parallel computation
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  LOG:        20111102 created
##
##  FUNCTIONS:  GLMM_MCMCwrapper
##
## ======================================================================

## *************************************************************
## GLMM_MCMCwrapper
## *************************************************************
GLMM_MCMCwrapper <- function(chain=1, data,
                             prior.alpha, init.alpha,
                             scale.b, prior.b, init.b,
                             prior.eps, init.eps,
                             Cpar, nMCMC, store, keep.chains, silent)
{
  thispackage <- "mixAK"

  
  ########## ========== Parameters from inits ========== ##########
  ########## =========================================== ##########
  Csigma_eps    <- init.eps[[chain]]$sigma
  CgammaInv_eps <- init.eps[[chain]]$gammaInv
  if (data$dimb){  
    CK_b <- init.b[[chain]]$K
    Cw_b <- c(init.b[[chain]]$w, rep(0, prior.b$Kmax - init.b[[chain]]$K))
    if (data$dimb == 1){
      Cmu_b <- c(init.b[[chain]]$mu, rep(0, prior.b$Kmax - init.b[[chain]]$K))
      CLi_b <- c(init.b[[chain]]$Li, rep(0, prior.b$Kmax - init.b[[chain]]$K))
    
    }else{
      Cmu_b <- c(t(init.b[[chain]]$mu), rep(0, data$dimb*(prior.b$Kmax - init.b[[chain]]$K)))
      CLi_b <- c(init.b[[chain]]$Li, rep(0, data$LTb*(prior.b$Kmax - init.b[[chain]]$K)))
    }
    CgammaInv_b <- init.b[[chain]]$gammaInv
    Cdf_b <- init.b[[chain]]$df
    Cr_b  <- init.b[[chain]]$r - 1    
    Cbb   <- as.numeric(t(init.b[[chain]]$b))
  }else{
    CK_b        <- 0
    Cw_b        <- 0
    Cmu_b       <- 0
    CLi_b       <- 0
    CgammaInv_b <- 0
    Cdf_b       <- 0
    Cr_b        <- 0
    Cbb         <- 0
  }  
  Calpha <- init.alpha[[chain]]
    
  
  ########## ========== Some additional parameters ##########
  ########## ===================================== ##########
  if (prior.b$priorK == "fixed") lsum_Ir_b <- Cpar$I * CK_b
  else                           lsum_Ir_b <- 1

  CshiftScale_b <- c(scale.b$shift, scale.b$scale)

  
########## ========== MCMC simulation                              ========== ##########
########## ================================================================== ##########
  if (!silent){
    cat(paste("\nChain number ", chain, "\n==============\n", sep=""))  
    cat(paste("MCMC sampling started on ", date(), ".\n", sep=""))
  }

  MCMC <- .C("GLMM_MCMC",
             Y_c                       = as.double(Cpar$Y_c),
             Y_d                       = as.integer(Cpar$Y_d),
             nonSilent_keepChain_nMCMC_R_cd_dist = as.integer(c(as.integer(!silent), store, nMCMC, Cpar$R_cd, Cpar$dist)),
             I_n                       = as.integer(c(Cpar$I, Cpar$n)),
             X                         = as.double(Cpar$X),
             #XtX                       = as.double(ifit$CXtX),               ### REMOVED ON 21/10/2009, XtX is computed directly in C++
             Z                         = as.double(Cpar$Z),
             #ZitZi                     = as.double(ifit$CZitZi),             ### REMOVED ON 20/10/2009, ZitZi is computed directly in C++
             p_fI_q_rI                 = as.integer(Cpar$p_fI_q_rI),
             shiftScale_b              = as.double(CshiftScale_b),
             priorDouble_eps           = as.double(Cpar$priorDouble_eps),
             priorInt_b                = as.integer(Cpar$priorInt_b),
             priorDouble_b             = as.double(Cpar$priorDouble_b),
             priorDouble_alpha         = as.double(Cpar$priorDouble_alpha),
             tune_scale_alpha_b        = as.double(Cpar$tune_scale_alpha_b),
             sigma_eps                 = as.double(Csigma_eps),
             gammaInv_eps              = as.double(CgammaInv_eps),
             K_b                       = as.integer(CK_b),
             w_b                       = as.double(Cw_b),
             mu_b                      = as.double(Cmu_b),
             Q_b                       = double(ifelse(data$dimb, data$LTb * prior.b$Kmax, 1)),
             Sigma_b                   = double(ifelse(data$dimb, data$LTb * prior.b$Kmax, 1)),
             Li_b                      = as.double(CLi_b),
             gammaInv_b                = as.double(CgammaInv_b),
             df_b                      = as.double(Cdf_b),
             r_b                       = as.integer(Cr_b),
             r_b_first                 = integer(Cpar$I),
             alpha                      = as.double(Calpha),
             b                         = as.double(Cbb),
             b_first                   = double(length(Cbb)),
             chsigma_eps               = double(ifelse(Cpar$R_cd["R_c"], Cpar$R_cd["R_c"] * nMCMC["keep"], 1)),
             chgammaInv_eps            = double(ifelse(Cpar$R_cd["R_c"], Cpar$R_cd["R_c"] * nMCMC["keep"], 1)),
             chK_b                     = integer(ifelse(data$dimb, nMCMC["keep"], 1)),
             chw_b                     = double(ifelse(data$dimb, prior.b$Kmax * nMCMC["keep"], 1)),
             chmu_b                    = double(ifelse(data$dimb, data$dimb * prior.b$Kmax * nMCMC["keep"], 1)),
             chQ_b                     = double(ifelse(data$dimb, data$LTb * prior.b$Kmax * nMCMC["keep"], 1)),
             chSigma_b                 = double(ifelse(data$dimb, data$LTb * prior.b$Kmax * nMCMC["keep"], 1)),
             chLi_b                    = double(ifelse(data$dimb, data$LTb * prior.b$Kmax * nMCMC["keep"], 1)),
             chgammaInv_b              = double(ifelse(data$dimb, data$dimb * nMCMC["keep"], 1)),
             chdf_b                    = double(ifelse(data$dimb, prior.b$Kmax * nMCMC["keep"], 1)),             
             chorder_b                 = integer(ifelse(data$dimb, prior.b$Kmax * nMCMC["keep"], 1)),
             chrank_b                  = integer(ifelse(data$dimb, prior.b$Kmax * nMCMC["keep"], 1)),
             chMeanData_b              = double(ifelse(data$dimb, data$dimb * nMCMC["keep"], 1)),
             chCorrData_b              = double(ifelse(data$dimb, data$LTb * nMCMC["keep"], 1)),
             chalpha                    = double(ifelse(data$lalpha, data$lalpha * nMCMC["keep"], 1)),
             chb                       = double(ifelse(data$dimb, ifelse(store["b"], Cpar$I * data$dimb * nMCMC["keep"], Cpar$I * data$dimb), 1)),
             chGLMMLogL                = double(nMCMC["keep"]),
             chLogL                    = double(nMCMC["keep"]),
             naccept_alpha              = integer(Cpar$R_cd["R_c"] + Cpar$R_cd["R_d"]),
             naccept_b                 = integer(Cpar$I),
             pm_eta_fixed              = double(Cpar$sumCn),
             pm_eta_random             = double(Cpar$sumCn),
             pm_meanY                  = double(Cpar$sumCn),
             pm_stres                  = double(Cpar$sumCn),
             pm_b                      = double(ifelse(data$dimb, data$dimb * Cpar$I, 1)),
             pm_w_b                    = double(ifelse(data$dimb, prior.b$Kmax, 1)),
             pm_mu_b                   = double(ifelse(data$dimb, data$dimb * prior.b$Kmax, 1)),
             pm_Q_b                    = double(ifelse(data$dimb, data$LTb * prior.b$Kmax, 1)),
             pm_Sigma_b                = double(ifelse(data$dimb, data$LTb * prior.b$Kmax, 1)),
             pm_Li_b                   = double(ifelse(data$dimb, data$LTb * prior.b$Kmax, 1)),
             pm_indGLMMLogL            = double(Cpar$I),
             pm_indLogL                = double(Cpar$I),
             pm_indLogpb               = double(Cpar$I),
             sum_Ir_b                  = integer(lsum_Ir_b),
             sum_Pr_b_b                = double(lsum_Ir_b),
             iter                      = as.integer(0),
             err                       = as.integer(0),
             PACKAGE=thispackage)            
  if (!silent) cat(paste("MCMC sampling finished on ", date(), ".\n", sep=""))
  if (MCMC$err) stop("Something went wrong.")


  ########## ========== State of MCMC (last and first kept) ========== ##########
  ########## ========================================================= ##########
  if (data$dimb){
    state.w_b       <- as.numeric(MCMC$w_b[1:MCMC$K_b])
    state_first.w_b <- as.numeric(MCMC$chw_b[1:MCMC$chK_b[1]])    
    names(state.w_b)       <- paste("w", 1:MCMC$K_b, sep="")
    names(state_first.w_b) <- paste("w", 1:MCMC$chK_b[1], sep="")    
    
    state.r_b       <- as.numeric(MCMC$r_b + 1)
    state_first.r_b <- as.numeric(MCMC$r_b_first + 1)    
    names(state.r_b) <- names(state_first.r_b) <- paste("r", 1:Cpar$I, sep="")
    
    state.gammaInv_b       <- as.numeric(MCMC$gammaInv_b)
    state_first.gammaInv_b <- as.numeric(MCMC$chgammaInv_b[1:data$dimb])    
    names(state_first.gammaInv_b) <- names(state.gammaInv_b) <- paste("gammaInv", 1:data$dimb, sep="")    
    
    if (data$dimb == 1){
      state.mu_b       <- as.numeric(MCMC$mu_b[1:MCMC$K_b])
      state_first.mu_b <- as.numeric(MCMC$chmu_b[1:MCMC$chK_b[1]])      
      names(state.mu_b)       <- paste("mu", 1:MCMC$K_b, sep="")
      names(state_first.mu_b) <- paste("mu", 1:MCMC$chK_b[1], sep="")      
      
      state.Li_b       <- as.numeric(MCMC$Li_b[1:MCMC$K_b])
      state_first.Li_b <- as.numeric(MCMC$chLi_b[1:MCMC$chK_b[1]])      
      names(state.Li_b)       <- paste("Li", 1:MCMC$K_b, sep="")
      names(state_first.Li_b) <- paste("Li", 1:MCMC$chK_b[1], sep="")      
      
      state.Sigma_b       <- (1 / state.Li_b)^2
      state_first.Sigma_b <- (1 / state_first.Li_b)^2      
      names(state.Sigma_b)       <- paste("Sigma", 1:MCMC$K_b, sep="")
      names(state_first.Sigma_b) <- paste("Sigma", 1:MCMC$chK_b[1], sep="")      
      
      state.Q_b       <- as.numeric(MCMC$Q_b[1:MCMC$K_b])
      state_first.Q_b <- as.numeric(MCMC$chQ_b[1:MCMC$chK_b[1]])      
      names(state.Q_b)       <- paste("Q", 1:MCMC$K_b, sep="")
      names(state_first.Q_b) <- paste("Q", 1:MCMC$chK_b[1], sep="")      

      state.b       <- as.numeric(MCMC$b)
      state_first.b <- as.numeric(MCMC$b_first)      
      names(state.b) <- names(state_first.b) <- 1:Cpar$I
    }else{
      state.mu_b       <- matrix(MCMC$mu_b[1:(data$dimb*MCMC$K_b)], ncol=data$dimb, byrow=TRUE)
      state_first.mu_b <- matrix(MCMC$chmu_b[1:(data$dimb*MCMC$chK_b[1])], ncol=data$dimb, byrow=TRUE)      
      rownames(state.mu_b)       <- paste("j", 1:MCMC$K_b, sep="")
      rownames(state_first.mu_b) <- paste("j", 1:MCMC$chK_b[1], sep="")      
      colnames(state.mu_b) <- colnames(state_first.mu_b) <- paste("m", 1:data$dimb, sep="")
      
      state.Li_b       <- as.numeric(MCMC$Li_b[1:(data$LTb*MCMC$K_b)])
      state_first.Li_b <- as.numeric(MCMC$chLi_b[1:(data$LTb*MCMC$chK_b[1])])      
      names(state.Li_b)       <- paste("Li", rep(1:MCMC$K_b, each=data$LTb), rep(data$naamLTb, MCMC$K_b), sep="")
      names(state_first.Li_b) <- paste("Li", rep(1:MCMC$chK_b[1], each=data$LTb), rep(data$naamLTb, MCMC$chK_b[1]), sep="")      
      
      state.Sigma_b       <- matrix(NA, ncol=data$dimb, nrow=data$dimb*MCMC$K_b)
      rownames(state.Sigma_b) <- paste("j", rep(1:MCMC$K_b, each=data$dimb), ".", rep(1:data$dimb, MCMC$K_b), sep="")
      colnames(state.Sigma_b) <- paste("m", 1:data$dimb, sep="")      
      for (j in 1:MCMC$K_b){
        tmpSigma <- matrix(0, nrow=data$dimb, ncol=data$dimb)
        tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- state.Li_b[((j-1)*data$LTb+1):(j*data$LTb)]
        tmpSigma <- tmpSigma %*% t(tmpSigma)
        tmpSigma <- chol2inv(chol(tmpSigma))
        state.Sigma_b[((j-1)*data$dimb+1):(j*data$dimb),] <- tmpSigma
      }

      state_first.Sigma_b <- matrix(NA, ncol=data$dimb, nrow=data$dimb*MCMC$chK_b[1])      
      rownames(state_first.Sigma_b) <- paste("j", rep(1:MCMC$chK_b[1], each=data$dimb), ".", rep(1:data$dimb, MCMC$chK_b[1]), sep="")      
      colnames(state_first.Sigma_b) <- paste("m", 1:data$dimb, sep="")      
      for (j in 1:MCMC$chK_b[1]){
        tmpSigma <- matrix(0, nrow=data$dimb, ncol=data$dimb)
        tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- state_first.Li_b[((j-1)*data$LTb+1):(j*data$LTb)]
        tmpSigma <- tmpSigma %*% t(tmpSigma)
        tmpSigma <- chol2inv(chol(tmpSigma))
        state_first.Sigma_b[((j-1)*data$dimb+1):(j*data$dimb),] <- tmpSigma
      }
      
      state.Q_b       <- as.numeric(MCMC$Q_b[1:(data$LTb*MCMC$K_b)])
      state_first.Q_b <- as.numeric(MCMC$chQ_b[1:(data$LTb*MCMC$chK_b[1])])      
      names(state.Q_b)       <- paste("Q", rep(1:MCMC$K_b, each=data$LTb), rep(data$naamLTb, MCMC$K_b), sep="")
      names(state_first.Q_b) <- paste("Q", rep(1:MCMC$chK_b[1], each=data$LTb), rep(data$naamLTb, MCMC$chK_b[1]), sep="")
      
      state.b       <- matrix(MCMC$b, ncol=data$dimb, nrow=Cpar$I, byrow=TRUE)
      state_first.b <- matrix(MCMC$b_first, ncol=data$dimb, nrow=Cpar$I, byrow=TRUE)      
      colnames(state.b) <- colnames(state_first.b) <- paste("b", 1:data$dimb, sep="")
      rownames(state.b) <- rownames(state_first.b) <- 1:Cpar$I
      
    }
    nCompTotal_b<- sum(MCMC$chK_b)
    freqK_b <- table(MCMC$chK_b)
    propK_b <- prop.table(freqK_b)
  }else{
    state.w_b <- state.r_b <- state.gamma_b <- state.mu_b <- state.Li_b <- state.Sigma_b <- state.Q_b <- state.b <- 0
    state_first.w_b <- state_first.r_b <- state_first.gamma_b <- state_first.mu_b <- state_first.Li_b <- state_first.Sigma_b <- state_first.Q_b <- state_first.b <- 0    
  }  

  if (data$lalpha){
    state.alpha       <- as.numeric(MCMC$alpha)
    state_first.alpha <- as.numeric(MCMC$chalpha[1:data$lalpha])    
    names(state.alpha) <- names(state_first.alpha) <- paste("alpha", 1:data$lalpha, sep="")
  }else{
    state.alpha <- state_first.alpha <- 0
  }  

  if (Cpar$R_cd["R_c"]){
    state.sigma_eps       <- as.numeric(MCMC$sigma_eps)
    state_first.sigma_eps <- as.numeric(MCMC$chsigma_eps[1:Cpar$R_cd["R_c"]])    
    names(state.sigma_eps) <- names(state_first.sigma_eps) <- paste("sigma", 1:Cpar$R_cd["R_c"], sep="")
    
    state.gammaInv_eps <- as.numeric(MCMC$gammaInv_eps)
    state_first.gammaInv_eps <- as.numeric(MCMC$chgammaInv_eps[1:Cpar$R_cd["R_c"]])    
    names(state.gammaInv_eps) <- names(state_first.gammaInv_eps) <- paste("gammaInv", 1:Cpar$R_cd["R_c"], sep="")
  }else{
    state.sigma_eps <- state.gammaInv_eps <- 0
    state_first.sigma_eps <- state_first.gammaInv_eps <- 0    
  }  

  
  ########## ========== Performance of MCMC ========== ##########
  ########## ========================================= ##########
  prop.accept.alpha <- MCMC$naccept_alpha / (nMCMC["keep"] * nMCMC["thin"])
  if (data$R > 1) names(prop.accept.alpha) <- data$name.response
  prop.accept.b <- MCMC$naccept_b / (nMCMC["keep"] * nMCMC["thin"])  
    
  
  ########## ========== Create a list to be returned ========== ##########
  ########## ================================================== ##########
  RET <- list(iter             = MCMC$iter,
              nMCMC            = nMCMC,
              dist             = data$dist,
              R                = c(Rc=as.numeric(Cpar$R_cd["R_c"]), Rd=as.numeric(Cpar$R_cd["R_d"])),
              p                = data$p,
              q                = data$q,
              fixed.intercept  = data$fixed.intercept,
              random.intercept = data$random.intercept,
              lalpha           = data$lalpha,
              dimb             = data$dimb,
              prior.alpha      = prior.alpha,
              prior.b          = prior.b,
              prior.eps        = prior.eps)


  if (data$lalpha){
    RET$init.alpha        <- init.alpha[[chain]]
    RET$state.first.alpha <- state_first.alpha    
    RET$state.last.alpha  <- state.alpha
    RET$prop.accept.alpha <- prop.accept.alpha    
  }  
  
  if (data$dimb){
    RET$init.b  <- init.b[[chain]]
    RET$state.first.b <- list(b        = state_first.b,               
                              K        = as.numeric(MCMC$chK_b[1]),  
                              w        = state_first.w_b,             
                              mu       = state_first.mu_b,            
                              Sigma    = state_first.Sigma_b,         
                              Li       = state_first.Li_b,            
                              Q        = state_first.Q_b,             
                              gammaInv = state_first.gammaInv_b,      
                              r        = state_first.r_b)    
    RET$state.last.b <- list(b        = state.b,               
                             K        = as.numeric(MCMC$K_b),  
                             w        = state.w_b,             
                             mu       = state.mu_b,            
                             Sigma    = state.Sigma_b,         
                             Li       = state.Li_b,            
                             Q        = state.Q_b,             
                             gammaInv = state.gammaInv_b,      
                             r        = state.r_b)
    RET$prop.accept.b <- prop.accept.b
    RET$scale.b <- scale.b
    RET$freqK_b <- freqK_b
    RET$propK_b <- propK_b    
  }                                 

  if (Cpar$R_cd["R_c"]){
    RET$init.eps <- init.eps[[chain]]
    RET$state.first.eps <- list(sigma    = state_first.sigma_eps,
                                gammaInv = state_first.gammaInv_eps)
    RET$state.last.eps <- list(sigma    = state.sigma_eps,
                               gammaInv = state.gammaInv_eps)    
  }                                       
  
  ########## ========== Posterior means of quantities computed in C++ ========== ##########
  ########## =================================================================== ##########
  RET$poster.mean.y <- list()
  used <- 0
  s    <- 1
  while (s <= Cpar$R_cd["R_c"]){
    ns    <- Cpar$n[((s-1)*Cpar$I+1):(s*Cpar$I)]
    index <- (used+1):(used + sum(ns))    
    used  <- index[length(index)]
    RET$poster.mean.y[[s]] <- data.frame(id         = rep(1:Cpar$I, ns),
                                         observed   = Cpar$Y_c[index],
                                         fitted     = as.numeric(MCMC$pm_meanY[index]),
                                         stres      = as.numeric(MCMC$pm_stres[index]),
                                         eta.fixed  = as.numeric(MCMC$pm_eta_fixed[index]),
                                         eta.random = as.numeric(MCMC$pm_eta_random[index]))
    s <- s + 1
  }
  used2 <- 0
  while (s <= Cpar$R_cd["R_c"] + Cpar$R_cd["R_d"]){
    ns     <- Cpar$n[((s-1)*Cpar$I+1):(s*Cpar$I)]
    index  <- (used+1):(used + sum(ns))    
    used   <- index[length(index)]
    index2 <- (used2+1):(used2 + sum(ns))    
    used2  <- index2[length(index2)]  
    RET$poster.mean.y[[s]] <- data.frame(id         = rep(1:Cpar$I, ns),
                                         observed   = Cpar$Y_d[index2],
                                         fitted     = as.numeric(MCMC$pm_meanY[index]),
                                         stres      = as.numeric(MCMC$pm_stres[index]),
                                         eta.fixed  = as.numeric(MCMC$pm_eta_fixed[index]),
                                         eta.random = as.numeric(MCMC$pm_eta_random[index]))
    s <- s + 1
  }  
  names(RET$poster.mean.y) <- colnames(data$y)
  
  if (data$dimb){
    MCMC$pm_b <- matrix(MCMC$pm_b, ncol=data$dimb, byrow=TRUE)
    RET$poster.mean.profile <- as.data.frame(MCMC$pm_b)
    colnames(RET$poster.mean.profile) <- paste("b", 1:data$dimb, sep="")

    RET$poster.mean.profile$Logpb         <- as.numeric(MCMC$pm_indLogpb)    
    RET$poster.mean.profile$Cond.Deviance <- as.numeric(-2 * MCMC$pm_indLogL)
    RET$poster.mean.profile$Deviance      <- as.numeric(-2 * MCMC$pm_indGLMMLogL)
    
    if (prior.b$priorK == "fixed"){
                                                      ##### I am not sure whether the posterior means (especially of variance components) are useful!
                                                      ##### In any case, they should be used with care
                                                      ##### -----------------------------------------------------------------------------------------
      RET$poster.mean.w_b <- as.numeric(MCMC$pm_w_b)
      names(RET$poster.mean.w_b) <- paste("w", 1:prior.b$Kmax, sep="")

      RET$poster.mean.mu_b <- matrix(MCMC$pm_mu_b, nrow=prior.b$Kmax, ncol=data$dimb, byrow=TRUE)
      rownames(RET$poster.mean.mu_b) <- paste("j", 1:prior.b$Kmax, sep="")
      colnames(RET$poster.mean.mu_b) <- paste("m", 1:data$dimb, sep="")

      RET$poster.mean.Q_b <- RET$poster.mean.Sigma_b <- RET$poster.mean.Li_b <- list()
      for (j in 1:prior.b$Kmax){
        tmpQ <- matrix(0, nrow=data$dimb, ncol=data$dimb)
        tmpQ[lower.tri(tmpQ, diag=TRUE)] <- MCMC$pm_Q_b[((j-1)*data$LTb+1):(j*data$LTb)]
        tmpQ[upper.tri(tmpQ, diag=FALSE)] <- t(tmpQ)[upper.tri(t(tmpQ), diag=FALSE)]
        RET$poster.mean.Q_b[[j]] <- tmpQ
      
        tmpSigma <- matrix(0, nrow=data$dimb, ncol=data$dimb)
        tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- MCMC$pm_Sigma_b[((j-1)*data$LTb+1):(j*data$LTb)]
        tmpSigma[upper.tri(tmpSigma, diag=FALSE)] <- t(tmpSigma)[upper.tri(t(tmpSigma), diag=FALSE)]
        RET$poster.mean.Sigma_b[[j]] <- tmpSigma
      
        tmpLi <- matrix(0, nrow=data$dimb, ncol=data$dimb)
        tmpLi[lower.tri(tmpLi, diag=TRUE)] <- MCMC$pm_Li_b[((j-1)*data$LTb+1):(j*data$LTb)]
        RET$poster.mean.Li_b[[j]] <- tmpLi      
      }
      names(RET$poster.mean.Q_b) <- names(RET$poster.mean.Sigma_b) <- names(RET$poster.mean.Li_b) <- paste("j", 1:prior.b$Kmax, sep="")    
    }      
  }else{
    RET$poster.mean.profile <- data.frame(LogL     = as.numeric(MCMC$pm_indLogL),
                                          Deviance = as.numeric(-2 * MCMC$pm_indGLMMLogL))
  }  

  
  ########## ========== Clustering based on posterior P(alloc = k | y) or on P(alloc = k | theta, b, y)    ========== ##########
  ########## ======================================================================================================== ##########
  if (data$dimb){
    if (prior.b$priorK == "fixed"){
      if (CK_b == 1){
        RET$poster.comp.prob_u <- RET$poster.comp.prob_b <- matrix(1, nrow = Cpar$I, ncol = 1)
      }else{

        ### Using mean(I(r=k))
        MCMC$sum_Ir_b <- matrix(MCMC$sum_Ir_b, ncol = CK_b, nrow = Cpar$I, byrow = TRUE)
        Denom <- apply(MCMC$sum_Ir_b, 1, sum)
        RET$poster.comp.prob_u <- MCMC$sum_Ir_b / matrix(rep(Denom, CK_b), ncol = CK_b, nrow = Cpar$I)

        ### Using mean(P(r=k | theta, b, y))
        MCMC$sum_Pr_b_b<- matrix(MCMC$sum_Pr_b_b, ncol = CK_b, nrow = Cpar$I, byrow = TRUE)
        RET$poster.comp.prob_b <- MCMC$sum_Pr_b_b/ matrix(rep(Denom, CK_b), ncol = CK_b, nrow = Cpar$I)        
      }  
    }  
  }  
  
  
  ########## ========== Additional posterior summaries                ========== ##########
  ########## =================================================================== ##########
  qProbs <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
  nSumm <- c("Mean", "Std.Dev.", "Min.", "2.5%", "1st Qu.", "Median", "3rd Qu.", "97.5%", "Max.")

  mean.Deviance  <- -2 * mean(MCMC$chGLMMLogL, na.rm=TRUE)
  quant.Deviance <- 2 * quantile(-MCMC$chGLMMLogL, prob=qProbs, na.rm=TRUE)  
  sd.Deviance    <- 2 * sd(MCMC$chGLMMLogL, na.rm=TRUE)
  summ.Deviance  <-  c(mean.Deviance, sd.Deviance, quant.Deviance)
  mean.Cond.Deviance  <- -2 * mean(MCMC$chLogL, na.rm=TRUE)
  quant.Cond.Deviance <- 2 * quantile(-MCMC$chLogL, prob=qProbs, na.rm=TRUE)  
  sd.Cond.Deviance    <- 2 * sd(MCMC$chLogL, na.rm=TRUE)
  summ.Cond.Deviance  <-  c(mean.Cond.Deviance, sd.Cond.Deviance, quant.Cond.Deviance)  
  RET$summ.Deviance <- data.frame(Deviance = summ.Deviance, Cond.Deviance = summ.Cond.Deviance)
  rownames(RET$summ.Deviance) <- nSumm
    
  if (data$lalpha){
    MCMC$chalpha <- matrix(MCMC$chalpha, ncol=data$lalpha, byrow=TRUE)
    colnames(MCMC$chalpha) <- paste("alpha", 1:data$lalpha, sep="")
    
    if (data$lalpha == 1){
      mean.alpha  <- mean(MCMC$chalpha, na.rm=TRUE)
      quant.alpha <- quantile(MCMC$chalpha, prob=qProbs, na.rm=TRUE)
      sd.alpha    <- sd(as.numeric(MCMC$chalpha), na.rm=TRUE)
      RET$summ.alpha     <- c(mean.alpha, sd.alpha, quant.alpha)
      names(RET$summ.alpha) <- nSumm      
    }else{
      mean.alpha  <- apply(MCMC$chalpha, 2, mean, na.rm=TRUE)
      quant.alpha <- apply(MCMC$chalpha, 2, quantile, prob=qProbs, na.rm=TRUE)
      sd.alpha    <- apply(MCMC$chalpha, 2, sd, na.rm=TRUE)
      RET$summ.alpha <- rbind(mean.alpha, sd.alpha, quant.alpha)
      RET$summ.alpha <- as.data.frame(RET$summ.alpha) 
      rownames(RET$summ.alpha) <- nSumm            
    }
  }  

  if (data$dimb){
    MCMC$chMeanData_b <- matrix(MCMC$chMeanData_b, ncol=data$dimb, byrow=TRUE)
    MCMC$chCorrData_b <- matrix(MCMC$chCorrData_b, ncol=data$LTb, byrow=TRUE)
    colnames(MCMC$chMeanData_b) <- paste("b.Mean.", 1:data$dimb, sep="")
    colnames(MCMC$chCorrData_b) <- paste("b.Corr", data$naamLTb, sep="")
    colnames(MCMC$chCorrData_b)[((0:(data$dimb-1))*(2*data$dimb - (0:(data$dimb-1)) + 1))/2 + 1] <- paste("b.SD.", 1:data$dimb, sep="")  
    
    if (data$dimb == 1){
      meanb.Mean  <- mean(MCMC$chMeanData_b, na.rm=TRUE)
      quantb.Mean <- quantile(MCMC$chMeanData_b, prob=qProbs, na.rm=TRUE)    
      sdb.Mean    <- sd(as.numeric(MCMC$chMeanData_b), na.rm=TRUE)
      RET$summ.b.Mean <- c(meanb.Mean, sdb.Mean, quantb.Mean)
      names(RET$summ.b.Mean) <- nSumm

      meanb.SDCorr  <- mean(MCMC$chCorrData_b, na.rm=TRUE)
      quantb.SDCorr <- quantile(MCMC$chCorrData_b, prob=qProbs, na.rm=TRUE)
      sdb.SDCorr    <- sd(as.numeric(MCMC$chCorrData_b), na.rm=TRUE)
      RET$summ.b.SDCorr <- c(meanb.SDCorr, sdb.SDCorr, quantb.SDCorr)
      names(RET$summ.b.SDCorr) <- nSumm    
    }else{
      meanb.Mean  <- apply(MCMC$chMeanData_b, 2, mean, na.rm=TRUE)
      quantb.Mean <- apply(MCMC$chMeanData_b, 2, quantile, prob=qProbs, na.rm=TRUE)    
      sdb.Mean    <- apply(MCMC$chMeanData_b, 2, sd, na.rm=TRUE)
      RET$summ.b.Mean <- rbind(meanb.Mean, sdb.Mean, quantb.Mean)
      rownames(RET$summ.b.Mean) <- nSumm

      meanb.SDCorr  <- apply(MCMC$chCorrData_b, 2, mean, na.rm=TRUE)
      quantb.SDCorr <- apply(MCMC$chCorrData_b, 2, quantile, prob=qProbs, na.rm=TRUE)    
      sdb.SDCorr    <- apply(MCMC$chCorrData_b, 2, sd, na.rm=TRUE)
      RET$summ.b.SDCorr <- rbind(meanb.SDCorr, sdb.SDCorr, quantb.SDCorr)
      rownames(RET$summ.b.SDCorr) <- nSumm        
    }  
  }  

  if (Cpar$R_cd["R_c"]){
    MCMC$chsigma_eps <- matrix(MCMC$chsigma_eps, ncol=Cpar$R_cd["R_c"], byrow=TRUE)
    colnames(MCMC$chsigma_eps) <- paste("sigma", 1:Cpar$R_cd["R_c"], sep="")

    if (Cpar$R_cd["R_c"] == 1){
      mean.sigma_eps  <- mean(MCMC$chsigma_eps, na.rm=TRUE)
      quant.sigma_eps <- quantile(MCMC$chsigma_eps, prob=qProbs, na.rm=TRUE)
      sd.sigma_eps    <- sd(as.numeric(MCMC$chsigma_eps), na.rm=TRUE)
      RET$summ.sigma_eps     <- c(mean.sigma_eps, sd.sigma_eps, quant.sigma_eps)
      names(RET$summ.sigma_eps) <- nSumm      
    }else{
      mean.sigma_eps  <- apply(MCMC$chsigma_eps, 2, mean, na.rm=TRUE)
      quant.sigma_eps <- apply(MCMC$chsigma_eps, 2, quantile, prob=qProbs, na.rm=TRUE)
      sd.sigma_eps    <- apply(MCMC$chsigma_eps, 2, sd, na.rm=TRUE)
      RET$summ.sigma_eps     <- rbind(mean.sigma_eps, sd.sigma_eps, quant.sigma_eps)
      rownames(RET$summ.sigma_eps) <- nSumm            
    }    
  }
  
  ########## ========== Chains for model parameters ========== ##########
  ########## ================================================= ##########
  if (keep.chains){
    RET$Deviance      <- as.numeric(-2 * MCMC$chGLMMLogL)
    RET$Cond.Deviance <- as.numeric(-2 * MCMC$chLogL)    
    
    if (data$dimb){
      ##### Chains for parameters of mixture distribution of b
      ##### -----------------------------------------------------
      RET$K_b <- as.numeric(MCMC$chK_b)
      MCMC$K_b <- NULL

      RET$w_b <- as.numeric(MCMC$chw_b[1:nCompTotal_b])
      MCMC$chw_b <- NULL

      RET$mu_b <- as.numeric(MCMC$chmu_b[1:(data$dimb*nCompTotal_b)])
      MCMC$chmu_b <- NULL

      RET$Li_b <- as.numeric(MCMC$chLi_b[1:(data$LTb*nCompTotal_b)])
      MCMC$chLi_b <- NULL    

      RET$Q_b <- as.numeric(MCMC$chQ_b[1:(data$LTb*nCompTotal_b)])
      MCMC$chQ_b <- NULL

      RET$Sigma_b <- as.numeric(MCMC$chSigma_b[1:(data$LTb*nCompTotal_b)])
      MCMC$chSigma_b <- NULL

      RET$gammaInv_b <- matrix(MCMC$chgammaInv_b, ncol=data$dimb, byrow=TRUE)
      colnames(RET$gammaInv_b) <- paste("gammaInv", 1:data$dimb, sep="")  
      MCMC$chgammaInv_b <- NULL
  
      RET$order_b <- as.numeric(MCMC$chorder_b[1:nCompTotal_b] + 1)
      MCMC$chorder_b <- NULL

      RET$rank_b <- as.numeric(MCMC$chrank_b[1:nCompTotal_b] + 1)
      MCMC$chrank_b <- NULL

      if (prior.b$priorK == "fixed"){
        RET$w_b <- matrix(RET$w_b, ncol=prior.b$Kmax, byrow=TRUE)
        colnames(RET$w_b) <- paste("w", 1:prior.b$Kmax, sep="")

        RET$mu_b <- matrix(RET$mu_b, ncol=data$dimb*prior.b$Kmax, byrow=TRUE)
        colnames(RET$mu_b) <- paste("mu.", rep(1:prior.b$Kmax, each=data$dimb), ".", rep(1:data$dimb, prior.b$Kmax), sep="")
    
        RET$Li_b <- matrix(RET$Li_b, ncol=data$LTb*prior.b$Kmax, byrow=TRUE)
        colnames(RET$Li_b) <- paste("Li", rep(1:prior.b$Kmax, each=data$LTb), rep(data$naamLTb, prior.b$Kmax), sep="")

        RET$Q_b <- matrix(RET$Q_b, ncol=data$LTb*prior.b$Kmax, byrow=TRUE)
        colnames(RET$Q_b) <- paste("Q", rep(1:prior.b$Kmax, each=data$LTb), rep(data$naamLTb, prior.b$Kmax), sep="")

        RET$Sigma_b <- matrix(RET$Sigma_b, ncol=data$LTb*prior.b$Kmax, byrow=TRUE)
        colnames(RET$Sigma_b) <- paste("Sigma", rep(1:prior.b$Kmax, each=data$LTb), rep(data$naamLTb, prior.b$Kmax), sep="")

        RET$order_b <- matrix(RET$order_b, ncol=prior.b$Kmax, byrow=TRUE)
        colnames(RET$order_b) <- paste("order", 1:prior.b$Kmax, sep="")

        RET$rank_b <- matrix(RET$rank_b, ncol=prior.b$Kmax, byrow=TRUE)
        colnames(RET$rank_b) <- paste("rank", 1:prior.b$Kmax, sep="")        
      }
            
      ##### Chains for characteristics of the mixture distribution of b
      ##### --------------------------------------------------------------
      RET$mixture_b <- as.data.frame(cbind(MCMC$chMeanData_b, MCMC$chCorrData_b))
      MCMC$chMeanData_b <- NULL
      MCMC$chCorrData_b <- NULL
      
      ##### Chains for random effects b
      ##### ------------------------------
      if (store["b"]){
        RET$b <- matrix(MCMC$chb, ncol=data$dimb*Cpar$I, byrow=TRUE)
        MCMC$chb <- NULL
        colnames(RET$b) <- paste("b.", rep(1:Cpar$I, each=data$dimb), ".", rep(1:data$dimb, Cpar$I), sep="")
      }        
    }

    if (data$lalpha){
      ##### Chains for regression coefficients alpha
      ##### ------------------------------------------
      RET$alpha <- MCMC$chalpha
      MCMC$chalpha <- NULL      
    }

    if (Cpar$R_cd["R_c"]){
      ##### Chains for parameters of distribution of residuals
      ##### -----------------------------------------------------
      RET$sigma_eps <- MCMC$chsigma_eps
      MCMC$chsigma_eps <- NULL

      RET$gammaInv_eps <- matrix(MCMC$chgammaInv_eps, ncol=Cpar$R_cd["R_c"], byrow=TRUE)
      colnames(RET$gammaInv_eps) <- paste("gammaInv", 1:Cpar$R_cd["R_c"], sep="")
      MCMC$chgammaInv_eps <- NULL
    }        
  }    

  
  ########## ========== Additional objects (added on 08/26/2010) ========== ##########
  ########## ============================================================== ##########
  RET$relabel_b <- list(type="mean", par=1)       #### default re-labeling is performed using the first margin of the mixture means
  RET$Cpar <- Cpar
  
  class(RET) <- "GLMM_MCMC"
  return(RET)      
}
