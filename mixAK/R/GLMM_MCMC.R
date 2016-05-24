##
##  PURPOSE:   Generalized linear mixed model with possibly several response variables
##             and normal mixtures in the distribution of random effects
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    06/07/2009
##              03/08/2009:  version for continuous responses working
##              10/11/2009:  version for combined discrete and continuous responses working
##              20/01/2011:  beta variables have been re-named to alpha variables
##              08/02/2013:  snow/snowfall support for parallel computation replaced by parallel package
##              01/06/2013:  minor changes following the review in JSS
##
##  FUNCTIONS:  GLMM_MCMC
##
## ================================================================================================

## *************************************************************
## GLMM_MCMC
## *************************************************************
##
GLMM_MCMC <- function(y, dist = "gaussian", id, x, z, random.intercept,
                      prior.alpha, init.alpha, init2.alpha,                      
                      scale.b,     prior.b,    init.b,      init2.b,
                      prior.eps,   init.eps,   init2.eps,
                      nMCMC = c(burn = 10, keep = 10, thin = 1, info = 10),
                      tuneMCMC = list(alpha = 1, b = 1),
                      store = c(b = FALSE), PED = TRUE, keep.chains = TRUE,
                      dens.zero = 1e-300, parallel = FALSE, cltype, silent = FALSE)
{
  #require("lme4")
  thispackage <- "mixAK"
  
  DEBUG <- FALSE

  silent <- as.logical(silent[1])
  if (is.na(silent)) silent <- FALSE

  EMin <- -100         ## exp(-(D1+D2)) = exp(-EMin) when computing importance sampling weights
                       ## which are equal to exp(-(D1+D2))
                       ## if D1 + D2 < EMin, where D1 = log(f(y|theta1)), D2 = log(f(y|theta2))
    ## -> these constants are passed to .C("GLMM_PED")
  
  
########## ========== Data ========== ##########
########## ========================== ##########
  dd <- GLMM_MCMCdata(y = y, dist = dist, id = id, x = x, z = z, random.intercept = random.intercept)
  rm(list=c("y", "dist", "id", "x", "z", "random.intercept"))
     ### use dd$y, dd$dist, dd$id, dd$x, dd$z, dd$random.intercept instead
     ### REMARK:  dd$x, dd$z are still without intercept column
  
  
########## ========== Initial fits ======================================== ##########
########## ========== and design information to be passed to C++ ========== ##########
########## ================================================================ ##########
  ifit <- GLMM_MCMCifit(do.init=TRUE, na.complete=FALSE,
                        y=dd$y, dist=dd$dist, id=dd$id, time=dd$time, x=dd$x, z=dd$z, random.intercept=dd$random.intercept,
                        xempty=dd$xempty, zempty=dd$zempty, Rc=dd$Rc, Rd=dd$Rd,
                        p=dd$p, p_fi=dd$p_fi, q=dd$q, q_ri=dd$q_ri, lalpha=dd$lalpha, dimb=dd$dimb)
  dd$x <- NULL
  dd$z <- NULL
     ### use ifit$x, ifit$z instead
     ### REMARK:  ifit$x, ifit$z contain intercept columns as well


########## ========== Prior distribution for fixed effects (alpha) ========== ##########
########## ================================================================= ##########
  palpha <- GLMM_MCMCprior.alpha(prior.alpha=prior.alpha, lalpha=dd$lalpha)
  rm(list="prior.alpha")                       ## use palpha$prior.alpha instead


########## ========== Prior distribution for error terms of gaussian responses ========== ##########
########## ============================================================================== ##########
  peps <- GLMM_MCMCprior.eps(prior.eps=prior.eps, Rc=dd$Rc, isigma=ifit$isigma, is.sigma=ifit$is.sigma)
  rm(list="prior.eps")                         ## use peps$prior.eps instead

  
########## ========== Shift and scale for random effects ========== ##########
########## ======================================================== ##########
  scale.b <- GLMM_MCMCscale.b(scale.b=scale.b, dimb=dd$dimb, iEranefVec=ifit$iEranefVec, iSDranefVec=ifit$iSDranefVec)
  
  
########## ========== Prior distribution for random effects ========== ##########
########## =========================================================== ##########
  pbb <- GLMM_MCMCprior.b(prior.b=prior.b, scale.b=scale.b, dimb=dd$dimb, iEranefVec=ifit$iEranefVec, iSDranefVec=ifit$iSDranefVec)
  rm(list="prior.b")                           ## use pbb$prior.b instead
  
  
########## ========== Initial values for fixed effects (alpha) ============== ##########
########## ================================================================= ##########
  init.alpha  <- GLMM_MCMCinit.alpha(init.alpha=init.alpha,  lalpha=dd$lalpha, ialpha=ifit$ialpha, number="")
  init2.alpha <- GLMM_MCMCinit.alpha(init.alpha=init2.alpha, lalpha=dd$lalpha, ialpha=ifit$ialpha2, number=2)

  
########## ========== Initial values for parameters related to the distribution of error terms of gaussian responses ========== ##########
########## ==================================================================================================================== #########
  init.eps  <- GLMM_MCMCinit.eps(init.eps=init.eps,  prior.eps=peps$prior.eps, Rc=dd$Rc, isigma=ifit$isigma[ifit$is.sigma],                        number="")
  init2.eps <- GLMM_MCMCinit.eps(init.eps=init2.eps, prior.eps=peps$prior.eps, Rc=dd$Rc, isigma=runif(dd$Rc, 0.9, 1.1)*ifit$isigma[ifit$is.sigma], number=2)  
  
  
########## ========== Initial values for parameters related to the distribution of random effects ========== ##########
########## ================================================================================================= ##########
  init.b <- GLMM_MCMCinit.b(init.b=init.b, prior.b=pbb$prior.b, scale.b=scale.b, 
                            id=dd$id, dimb=dd$dimb, LTb=dd$LTb, naamLTb=dd$naamLTb,
                            I=ifit$I, ibMat=ifit$ibMat, iEranefVec=ifit$iEranefVec, iSDranefVec=ifit$iSDranefVec, number="")
  init2.b <- GLMM_MCMCinit.b(init.b=init2.b, prior.b=pbb$prior.b, scale.b=scale.b, 
                             id=dd$id, dimb=dd$dimb, LTb=dd$LTb, naamLTb=dd$naamLTb,
                             I=ifit$I, ibMat=ifit$ibMat2,
                             iEranefVec=rnorm(length(ifit$iEranefVec), mean=ifit$iEranefVec, sd=1*ifit$iSEranefVec),
                             iSDranefVec=runif(length(ifit$iSDranefVec), 0.9, 1.1)*ifit$iSDranefVec,
                             number=2)
  
  
########## ========== nMCMC ========== ##########
########## =========================== ##########
  if (length(nMCMC) != 4) stop("nMCMC must be of length 4")
  if (is.null(names(nMCMC))) names(nMCMC) <- c("burn", "keep", "thin", "info")
  names.nMCMC <- names(nMCMC)
  if (!match("burn", names.nMCMC, nomatch=0)) stop(paste("nMCMC[", dQuote("burn"), "] must be specified", sep=""))
  else                                        n.burn <- nMCMC["burn"]
  if (!match("keep", names.nMCMC, nomatch=0)) stop(paste("nMCMC[", dQuote("keep"), "] must be specified", sep=""))
  else                                        n.keep <- nMCMC["keep"]
  if (!match("thin", names.nMCMC, nomatch=0)) stop(paste("nMCMC[", dQuote("thin"), "] must be specified", sep=""))
  else                                        n.thin <- nMCMC["thin"]
  if (!match("info", names.nMCMC, nomatch=0)) stop(paste("nMCMC[", dQuote("info"), "] must be specified", sep=""))
  else                                        n.info <- nMCMC["info"]
  nMCMC <- c(n.burn, n.keep, n.thin, n.info)
  names(nMCMC) <- c("burn", "keep", "thin", "info")  
  if (nMCMC["burn"] < 0) stop(paste("nMCMC[", dQuote("burn"), "] must be non-negative", sep=""))
  if (nMCMC["keep"] <= 0) stop(paste("nMCMC[", dQuote("keep"), "] must be positive", sep=""))  
  if (nMCMC["thin"] <= 0) stop(paste("nMCMC[", dQuote("thin"), "] must be positive", sep=""))
  if (nMCMC["info"] <= 0 | nMCMC["info"] > max(nMCMC["burn"], nMCMC["keep"])) nMCMC["info"] <- max(nMCMC["burn"], nMCMC["keep"])


########## ========== tuneMCMC ========== ##########
########## ============================== ##########
  if (!is.list(tuneMCMC)) stop("tuneMCMC must be a list")
  intuneMCMC <- names(tuneMCMC)
  itune.alpha <- match("alpha", intuneMCMC, nomatch=NA)
  itune.b    <- match("b", intuneMCMC, nomatch=NA)

  tuneMCMC_alpha_b <- numeric()
  if (dd$Rd){
    if (is.na(itune.alpha)) tuneMCMC$alpha <- rep(1, dd$Rd)
    if (length(tuneMCMC$alpha) == 1) tuneMCMC$alpha <- rep(tuneMCMC$alpha, dd$Rd)
    if (length(tuneMCMC$alpha) != dd$Rd) stop(paste("tuneMCMC$alpha must be of length ", dd$Rd, sep=""))
    if (any(is.na(tuneMCMC$alpha))) stop("NA in tuneMCMC$alpha")        
    if (any(tuneMCMC$alpha <= 0)) stop("tuneMCMC$alpha must be all positive")
    
    tuneMCMC_alpha_b <- c(tuneMCMC_alpha_b, tuneMCMC$alpha)
  }else{
    tuneMCMC$alpha <- 1
  }  
  
  if (dd$dimb){
    if (is.na(itune.b)) tuneMCMC$b <- 1
    if (length(tuneMCMC$b) != 1) stop(paste("tuneMCMC$b must be of length ", 1, sep=""))
    if (any(is.na(tuneMCMC$b))) stop("NA in tuneMCMC$b")        
    if (any(tuneMCMC$b <= 0)) stop("tuneMCMC$b must be all positive")

    tuneMCMC_alpha_b <- c(tuneMCMC_alpha_b, tuneMCMC$b)
  }else{    
    tuneMCMC$b <- 1
  }  

  if (length(tuneMCMC_alpha_b) == 0) tuneMCMC_alpha_b <- 1
  
  
########## ========== store ========== ##########
########## =========================== ##########
  if (length(store) != 1) stop("store must be of length 1")
  if (is.null(names(store))) names(store) <- c("b")
  names.store <- names(store)
  if (!match("b", names.store, nomatch=0)) stop(paste("store[", dQuote("b"), "] must be specified", sep=""))
  else                                     store.b <- store["b"]
  store <- c(store.b)
  names(store) <- c("b")
  if (!dd$dimb) store["b"] <- FALSE  


####### Parameters passed to C++ which will be stored also in the resulting object (to be able to use them in related functions)
####### ========================================================================================================================
  R_cd <- c(dd$Rc, dd$Rd)
  names(R_cd) <- c("R_c", "R_d")
  p_fI_q_rI <- c(dd$p, dd$CfixedIntcpt, dd$q, dd$CrandomIntcpt)
  names(p_fI_q_rI) <- paste(rep(c("p", "fixedIntcpt", "q", "randomIntcpt"), each=sum(R_cd)), rep(1:sum(R_cd), 4), sep="")
  
  Cpar <- list(Y_c                = ifit$Cy_c,
               Y_d                = ifit$Cy_d,
               R_cd               = R_cd,
               dist               = dd$ndist,
               I                  = ifit$I,
               n                  = ifit$Cn,
               sumCn              = ifit$sumCn,
               X                  = ifit$CX,
               Z                  = ifit$CZ,
               p_fI_q_rI          = p_fI_q_rI,
               priorDouble_eps    = peps$CpriorDouble_eps,
               priorInt_b         = pbb$CpriorInt_b,
               priorDouble_b      = pbb$CpriorDouble_b,
               priorDouble_alpha  = palpha$CpriorDouble_alpha,
               tune_scale_alpha   = tuneMCMC$alpha,
               tune_scale_b       = tuneMCMC$b,               
               tune_scale_alpha_b = tuneMCMC_alpha_b)
  ifit$Cy_c <- NULL
  ifit$Cy_d <- NULL
  ifit$C_n  <- NULL
  ifit$CX   <- NULL
  ifit$CZ   <- NULL  

  #return(list(Cpar=Cpar, scb=scb, nMCMC=nMCMC, store=store,
  #            Csigma_eps=Csigma_eps, CgammaInv_eps=CgammaInv_eps,
  #            CK_b=CK_b, Cw_b=Cw_b, Cmu_b=Cmu_b,
  #            dd=dd, prior.b=pbb$prior.b,
  #            CLi_b=CLi_b, CgammaInv_b=CgammaInv_b, Cr_b=Cr_b,
  #            ifit=ifit, Calpha=Calpha, Cbb=Cbb))

 
########## ========== Run MCMC ========== ##########
########## ============================== ##########  
  if (PED){
    if (parallel){
      #require("parallel")

      if (parallel::detectCores() < 2) warning("It does not seem that at least 2 CPU cores are available needed for efficient parallel generation of the two chains.")
      if (missing(cltype)) cl <- parallel::makeCluster(2) else cl <- parallel::makeCluster(2, type = cltype)
      if (!silent) cat(paste("Parallel MCMC sampling of two chains started on ", date(), ".\n", sep=""))      
      RET <- parallel::parLapply(cl, 1:2, GLMM_MCMCwrapper,
                                 data = dd,
                                 prior.alpha = palpha$prior.alpha, init.alpha = list(init.alpha, init2.alpha),
                                 scale.b = scale.b, prior.b = pbb$prior.b, init.b = list(init.b, init2.b),
                                 prior.eps = peps$prior.eps, init.eps = list(init.eps, init2.eps),
                                 Cpar = Cpar, nMCMC = nMCMC, store = store, keep.chains = keep.chains, silent = silent)
      if (!silent) cat(paste("Parallel MCMC sampling finished on ", date(), ".\n", sep=""))
      parallel::stopCluster(cl)
    }else{
      RET <- lapply(1:2, GLMM_MCMCwrapper,
                         data = dd,
                         prior.alpha = palpha$prior.alpha, init.alpha = list(init.alpha, init2.alpha),
                         scale.b = scale.b, prior.b = pbb$prior.b, init.b = list(init.b, init2.b),
                         prior.eps = peps$prior.eps, init.eps = list(init.eps, init2.eps),
                         Cpar = Cpar, nMCMC = nMCMC, store = store, keep.chains = keep.chains, silent = silent)
    }
    
    if (!silent) cat(paste("\nComputation of penalized expected deviance started on ", date(), ".\n", sep=""))
    resPED <- .C("GLMM_PED", PED                  = double(5),
                             pm.indDevObs         = double(Cpar$I),
                             pm.indpopt           = double(Cpar$I),
                             pm.windpopt          = double(Cpar$I),
                             invalid.indDevObs    = integer(Cpar$I),
                             invalid.indpopt      = integer(Cpar$I),
                             invalid.windpopt     = integer(Cpar$I),                   
                             sum.ISweight         = double(Cpar$I),
                             #ch.ISweight          = double(Cpar$I*nMCMC["keep"]),
                             chGLMMLogL1          = double(nMCMC["keep"]),
                             chGLMMLogL2          = double(nMCMC["keep"]),
                             chGLMMLogL_repl1_ch1 = double(nMCMC["keep"]),
                             chGLMMLogL_repl1_ch2 = double(nMCMC["keep"]),
                             chGLMMLogL_repl2_ch1 = double(nMCMC["keep"]),
                             chGLMMLogL_repl2_ch2 = double(nMCMC["keep"]),                 
                             err                  = integer(1),
                             Y_c                  = as.double(Cpar$Y_c),
                             Y_d                  = as.integer(Cpar$Y_d),
                             R_cd_dist            = as.integer(c(Cpar$R_cd, Cpar$dist)),
                             I_n                  = as.integer(c(Cpar$I, Cpar$n)),
                             X                    = as.double(Cpar$X),
                             Z                    = as.double(Cpar$Z),
                             p_fI_q_rI            = as.integer(Cpar$p_fI_q_rI),
                             distribution_b       = as.integer(Cpar$priorInt_b["distribution"]),
                             shiftScale_b         = if (dd$dimb) as.double(c(scale.b$shift, scale.b$scale)) else as.double(c(0, 1)),
                             chsigma_eps1         = if (dd$Rc) as.double(t(RET[[1]]$sigma_eps)) else as.double(0),
                             chK_b1               = if (dd$dimb) as.integer(RET[[1]]$K_b) else as.double(0),
                             chw_b1               = if (dd$dimb) as.double(t(RET[[1]]$w_b)) else as.double(0),
                             chmu_b1              = if (dd$dimb) as.double(t(RET[[1]]$mu_b)) else as.double(0),
                             chLi_b1              = if (dd$dimb) as.double(t(RET[[1]]$Li_b)) else as.double(0),
                             chQ_b1               = if (dd$dimb) as.double(t(RET[[1]]$Q_b)) else as.double(0),                 
                             chdf_b1              = if (dd$dimb & Cpar$priorInt_b["distribution"] == 1) as.double(t(RET[[1]]$df_b)) else as.double(0),
                             chbeta1              = if (dd$lalpha) as.double(t(RET[[1]]$alpha)) else as.double(0),
                             bhat1                = if (dd$dimb) as.double(t(RET[[1]]$poster.mean.profile[, 1:dd$dimb])) else as.double(0),
                             chsigma_eps2         = if (dd$Rc) as.double(t(RET[[2]]$sigma_eps)) else as.double(0),
                             chK_b2               = if (dd$dimb) as.integer(RET[[2]]$K_b) else as.double(0),
                             chw_b2               = if (dd$dimb) as.double(t(RET[[2]]$w_b)) else as.double(0),
                             chmu_b2              = if (dd$dimb) as.double(t(RET[[2]]$mu_b)) else as.double(0),
                             chLi_b2              = if (dd$dimb) as.double(t(RET[[2]]$Li_b)) else as.double(0),
                             chQ_b2               = if (dd$dimb) as.double(t(RET[[2]]$Q_b)) else as.double(0),                 
                             chdf_b2              = if (dd$dimb & Cpar$priorInt_b["distribution"] == 1) as.double(t(RET[[2]]$df_b)) else as.double(0),
                             chbeta2              = if (dd$lalpha) as.double(t(RET[[2]]$alpha)) else as.double(0),
                             bhat2                = if (dd$dimb) as.double(t(RET[[2]]$poster.mean.profile[, 1:dd$dimb])) else as.double(0),
                             M                    = as.integer(nMCMC["keep"]),
                             Dens_ZERO            = as.double(dens.zero),
                             EMin                 = as.double(EMin),
                 PACKAGE = thispackage)    
    if (!silent) cat(paste("Computation of penalized expected deviance finished on ", date(), ".\n\n", sep=""))
    if (resPED$err) stop("Something went wrong.")

    names(resPED$PED) <- c("D.expect", "p(opt)", "PED", "wp(opt)", "wPED")
    RET$PED <- resPED$PED    
    
    RET$D         <- resPED$pm.indDevObs
    RET$popt      <- resPED$pm.indpopt
    RET$wpopt     <- resPED$pm.windpopt
    RET$inv.D     <- resPED$invalid.indDevObs
    RET$inv.popt  <- resPED$invalid.indpopt
    RET$inv.wpopt <- resPED$invalid.windpopt    
    RET$sumISw    <- resPED$sum.ISweight
    #RET$ISw       <- matrix(resPED$ch.ISweight, nrow=dd$n, byrow=TRUE)

    RET$Deviance1 <- -2 * resPED$chGLMMLogL1
    RET$Deviance2 <- -2 * resPED$chGLMMLogL2
    RET$Deviance_repl1_ch1 <- -2 * resPED$chGLMMLogL_repl1_ch1
    RET$Deviance_repl1_ch2 <- -2 * resPED$chGLMMLogL_repl1_ch2
    RET$Deviance_repl2_ch1 <- -2 * resPED$chGLMMLogL_repl2_ch1
    RET$Deviance_repl2_ch2 <- -2 * resPED$chGLMMLogL_repl2_ch2            
    
    class(RET) <- "GLMM_MCMClist"    
  }else{    
    RET <- GLMM_MCMCwrapper(chain = 1, data = dd,
                            prior.alpha = palpha$prior.alpha, init.alpha = list(init.alpha),
                            scale.b = scale.b, prior.b = pbb$prior.b, init.b = list(init.b),
                            prior.eps = peps$prior.eps, init.eps = list(init.eps),
                            Cpar = Cpar, nMCMC = nMCMC, store = store, keep.chains = keep.chains, silent = silent)
  }
  
  return(RET)    
}

