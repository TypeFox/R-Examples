##
##  PURPOSE:   Longitudinal discriminant analysis
##             based on GLMM MCMC fits (with possibly several response variables)
##             * version 2
##             * Difference to GLMM_longitDA:
##                   - Only one prediction for each longitudinal profile is returned
##                     whereas with GLMM_longitDA, "sequential" predictions of a group
##                     membership are directly calculated.
##                   - GLMM_longitDA2 allows also for discrete responses
##                     which was not the case for GLMM_longitDA.
##                                      
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  LOG:       10/04/2015: start implementation allowing for discrete responses as well
##             
##  FUNCTIONS:  GLMM_longitDA2
##
## ================================================================================================

## *************************************************************
## GLMM_longitDA2
## *************************************************************
##
GLMM_longitDA2 <- function(mod, w.prior, y, id, x, z, xz.common = TRUE, keep.comp.prob = FALSE, level = 0.95, info, silent = FALSE)
{
  thispackage <- "mixAK" 

  
########## ========== Some input checks ========== ##########
########## ======================================= ##########
  ### This part is the same as in GLMM_longitDA on 10/04/2015.
  ###
  if (!is.list(mod)) stop("mod must be a list")
  if (length(mod) < 2) stop("mod must be a list of length 2 or more")
  if (any(sapply(mod, class) != "GLMM_MCMC")) stop("all components of mod must be of class GLMM_MCMC")
  nClust <- length(mod)

  if (length(w.prior) != nClust) stop(paste("w.prior must be of length ", nClust, sep=""))
  if (any(w.prior < 0)) stop("w.prior must all be non-negative")
  w.prior <- w.prior / sum(w.prior)
  
  Rc <- mod[[1]]$R["Rc"]
  Rd <- mod[[1]]$R["Rd"]
  R <- Rc + Rd

  keepMCMC <- mod[[1]]$nMCMC["keep"]
  for (cl in 2:nClust){
    if (mod[[cl]]$R["Rc"] != Rc) stop(paste("cluster number ", cl, " has ", mod[[cl]]$R["Rc"], " cont. responses (must be ", Rc, ")", sep=""))
    if (mod[[cl]]$R["Rd"] != Rd) stop(paste("cluster number ", cl, " has ", mod[[cl]]$R["Rd"], " discrete responses (must be ", Rd, ")", sep=""))
    keepMCMC <- c(keepMCMC, mod[[cl]]$nMCMC["keep"])
  }  

  
########## ========== Design matrices ============ ##########
########## ======================================= ##########
  ### This part is the same as in GLMM_longitDA on 10/04/2015.
  ###
  if (!xz.common){
    if (!is.list(x)) stop("x must be a list when xz.common is FALSE")
    if (length(x) != nClust) stop(paste("x must be a list of length ", nClust, sep=""))
    if (!is.list(z)) stop("z must be a list when xz.common is FALSE")
    if (length(z) != nClust) stop(paste("z must be a list of length ", nClust, sep=""))    
  }  


########## ========== Work out the data ============ ##########
########## ========================================= ##########
  ### This part is only more or less the same as
  ### in GLMM_longitDA on 10/04/2015. Some small
  ### parts differ here, however.
  ###  
  Cp <- Cq <- CfixedIntcpt <- CrandomIntcpt <- numeric()
  CX <- CZ <- numeric()
  CshiftScale_b <- numeric()
  Kmax_b <- numeric()
  chK_b <- chw_b <- chmu_b <- chLi_b <- numeric()
  chsigma_eps <- numeric()
  chalpha <- numeric()
  #
  #bnew.init <- numeric()
  #rbnew.init <- numeric()

  dimb_cl <- numeric(nClust)
  id_cl <- list()
  
  for (cl in 1:nClust){     ### Loop over clusters

    ##### Design stuff
    ##### ---------------
    if (xz.common){
      dd <- GLMM_MCMCdata(y = y, dist = mod[[cl]]$dist, id = id, x = x, z = z, random.intercept = mod[[cl]]$random.intercept)
    }else{
      dd <- GLMM_MCMCdata(y = y, dist = mod[[cl]]$dist, id = id, x = x[[cl]], z = z[[cl]], random.intercept = mod[[cl]]$random.intercept)
    }

    dimb_cl[cl] <- dd$dimb
    id_cl[[cl]] <- dd$id
    
    ifit <- GLMM_MCMCifit(do.init = FALSE, na.complete = FALSE,   ## in GLMM_longitDA, na.complete = TRUE
                          y = dd$y, dist = dd$dist, id = dd$id, x = dd$x, z = dd$z, random.intercept = dd$random.intercept,
                          time = dd$time, xempty = dd$xempty, zempty = dd$zempty, Rc = dd$Rc, Rd = dd$Rd,
                          p = dd$p, p_fi = dd$p_fi, q = dd$q, q_ri = dd$q_ri, lalpha = dd$lalpha, dimb = dd$dimb)
    ## ifit$x[[s]], ifit$z[[s]] now contains also intercept columns if these should be included
    ## * also rows corresponding to missing values are removed
    ## * there is the same number of observations for each cluster
    ## ifit$CX = non-intercept columns of X matrices,
    ##           for those s, where x = "empty", there is nothing in CX
    ##         = 0 if for all s, x = "empty"
    ## ifit$CZ = non-intercept columns of Z matrices,
    ##           for those s, where z = "empty", there is nothing in CZ
    ##         = 0 if for all s, z = "empty"
    
    ##### Initial values for random effects and component allocations
    ##### (within a cluster - model)
    ##### -------------------------------------------------------------
    #####
    ##### Finally not needed
    #####
    #lchcl <- length(mod[[cl]]$K)
    #iSigma.cl <- SP2Rect(mod[[cl]]$Sigma[lchcl, 1:dd$LTb], dd$dimb)
    #if (mod[[cl]]$K[lchcl] > 1){
    #  for (kcl in 2:mod[[cl]]$K[lchcl]){
    #    iSigma.cl <- rbind(iSigma.cl, SP2Rect(mod[[cl]]$Sigma[lchcl, ((kcl - 1) * dd$LTb + 1):(kcl * dd$LTb)], dd$dimb))
    #  }
    #}
    #init.b.cl <- list(K  = mod[[cl]]$K[lchcl],
    #                  w  = mod[[cl]]$w[lchcl,],
    #                  mu = matrix(mod[[cl]]$mu[lchcl,], ncol = dd$dimb, nrow = mod[[cl]]$K[lchcl], byrow = TRUE),
    #                  Sigma = iSigma.cl)    
    #if (dd$dimb == 1){
    #  iEranefVec.cl <- mod[[cl]]$summ.b.Mean["Mean"]
    #  iSDranefVec.cl <- mod[[cl]]$summ.b.SDCorr["Mean"]
    #  ibMat.cl <- matrix(rep(iEranefVec.cl, ifit$I), ncol = dd$dimb, byrow = TRUE)
    #}else{
    #  iEranefVec.cl <- mod[[cl]]$summ.b.Mean["Mean",]
    #  iSDranefVec.cl <- mod[[cl]]$summ.b.SDCorr["Mean", grep("b.SD", colnames(mod[[cl]]$summ.b.SDCorr))]
    #  ibMat.cl <- matrix(rep(iEranefVec.cl, ifit$I), ncol = dd$dimb, byrow = TRUE)
    #}    

      ### Initial values of random effects = posterior mean of the overall mixture mean.
      ### Also initial component allocations are computed here.
    #ibnew <- GLMM_MCMCinit.b(init.b = init.b.cl, prior.b = mod[[cl]]$prior.b, scale.b = mod[[cl]]$scale.b,
    #                         id = dd$id, dimb = dd$dimb, LTb = dd$LTb, naamLTb = dd$naamLTb,
    #                         I = ifit$I, ibMat = ibMat.cl, iEranefVec = iEranefVec.cl, iSDranefVec = iSDranefVec.cl, number = "")    
  
    #bnew.init  <- c(bnew.init, as.numeric(t(ibnew$b)))
    #rbnew.init <- c(rbnew.init, ibnew$r - 1)

    
    ##### Variables common to all clusters
    ##### ------------------------------------
    if (cl == 1){
      dist1 <- dd$dist
      Cdist <- dd$ndist

      Cy_c  <- ifit$Cy_c
      Cy_d  <- ifit$Cy_d
      I     <- ifit$I
      Cn    <- ifit$Cn
      sumCn <- ifit$sumCn

      ID <- ifit$ID[[1]]
      #TIME <- ifit$time

      nrow_pi <- sum(ifit$n[[1]])
    }  
    
    ##### Check whether supplied y is compatible with these used to fit mod
    ##### -------------------------------------------------------------------
    if (ncol(dd$y) != R) stop(paste("y must have ", R, " columns", sep=""))

    ##### Check whether common variables (across clusters) are really common
    ##### --------------------------------------------------------------------
    if (any(dd$dist != dist1)) stop(paste("dist of cluster ", cl, " is incompatible with that of cluster 1", sep=""))

    if (length(ifit$Cy_c) != length(Cy_c)) stop(paste("number of NA's for cluster ", cl, " leads to incompatibility in Cy_c", sep=""))
    if (length(ifit$Cy_d) != length(Cy_d)) stop(paste("number of NA's for cluster ", cl, " leads to incompatibility in Cy_d", sep=""))    
    if (ifit$I != I) stop(paste("I for cluster ", cl, " is incompatible with that of clulster 1", sep=""))
    if (length(ifit$Cn) != length(Cn)) stop(paste("Cn for cluster ", cl, " is incompatible with that of cluster 1", sep=""))
    if (any(ifit$Cn != Cn)) stop(paste("Cn for cluster ", cl, " is incompatible with that of cluster 1", sep=""))    
    
    ##### Check whether supplied x and z is compatible with these used to fit
    ##### ------------------------------------------------------------------------------
    for (s in 1:R){      
      if (dd$xempty[s]){
        if (mod[[cl]]$p[[s]] != 0) stop(paste("x[[", s, "]] should not be empty for cluster ", cl, sep=""))
      }else{
        if (ncol(dd$x[[s]]) != mod[[cl]]$p[s]) stop(paste("x[[", s, "]] matrix for cluster ", cl, " must have ", mod[[cl]]$p[s], " columns", sep=""))
      }

      if (dd$zempty[s]){
        if (dd$z[[s]] != "empty") stop(paste("z[[", s, "]] must be empty for cluster ", cl, sep=""))
      }else{
        if (ncol(dd$z[[s]]) != mod[[cl]]$q[s]) stop(paste("z[[", s, "]] matrix for cluster ", cl, " must have ", mod[[cl]]$q[s], " columns", sep=""))
      }

      ##if (any(ifit$n[[s]] != ifit$n[[1]])) stop("BUG in the function, contact AK")    ### This was required for GLMM_longitDA,
                                                                                        ### it is not required for GLMM_longitDA2.
    }
    
    ##### Vectors to be supplied to C++
    ##### ---------------------------------
    Cp <- c(Cp, dd$p)
    Cq <- c(Cq, dd$q)    
    CfixedIntcpt  <- c(CfixedIntcpt, dd$CfixedIntcpt)
    CrandomIntcpt <- c(CrandomIntcpt, dd$CrandomIntcpt)
    CX     <- c(CX, ifit$CX)
    CZ     <- c(CZ, ifit$CZ)

    CshiftScale_b <- c(CshiftScale_b, mod[[cl]]$scale.b$shift, mod[[cl]]$scale.b$scale)
    
    Kmax_b <- c(Kmax_b, mod[[cl]]$prior.b$Kmax)
    if (mod[[cl]]$prior.b$priorK == "fixed"){
      chK_b  <- c(chK_b, mod[[cl]]$K_b)
      chw_b  <- c(chw_b, t(mod[[cl]]$w_b))
      chmu_b <- c(chmu_b, t(mod[[cl]]$mu_b))    
      chLi_b <- c(chLi_b, t(mod[[cl]]$Li_b))
    }else{
      stop("not implemented for variable K")
    }  
    
    if (mod[[cl]]$lalpha)  chalpha <- c(chalpha, t(mod[[cl]]$alpha))
    if (mod[[cl]]$R["Rc"]) chsigma_eps <- c(chsigma_eps, t(mod[[cl]]$sigma_eps))
  }     ### End of loop over clusters

  if (!length(chalpha))     chalpha <- 0
  if (!length(chsigma_eps)) chsigma_eps <- 0  
  
  ### !!! CX, CZ contain one 0 when there is no such matrix
  ###     in a model for specific cluster.

  if (missing(info)) info <- min(keepMCMC)
  if (info <= 0) info <- 1

  Cdistribution_b <- 0     ### = mixture of normal distributions,
                           ### Nothing else than this is currently implemented...

  if (FALSE){
    cat("Cy_c:\n")
    print(Cy_c)

    cat("Cy_d:\n")
    print(Cy_d)

    cat("CX:\n")
    print(CX)

    cat("CZ:\n")
    print(CZ)

    cat("R_c = ", Rc, ", R_d = ", Rd, ", nClust = ", nClust, ", I = ", I, "\n", sep="")
    cat("Cdist:\n")
    print(Cdist)

    cat("Cn:\n")
    print(Cn)

    cat("Cp:\n")
    print(Cp)
    
    cat("CfixedIntcpt:\n")
    print(CfixedIntcpt)

    cat("Cq:\n")
    print(Cq)
    
    cat("CrandomIntcpt:\n")
    print(CrandomIntcpt)

    cat("CshiftScale_b:\n")
    print(CshiftScale_b)
        
    cat("keepMCMC:\n")
    print(keepMCMC)        
  }  
  
  
########## ========== Main calculation  ============ ##########
########## ========================================= ##########
  
    ##### Length of memory for logf_***_i
    ##### ---------------------------------
  if (keep.comp.prob){
    l_logf_i <- I * sum(keepMCMC)
  }else{
    l_logf_i <- I
  }    
  
    ##### Main calculation
    ##### ---------------------------------
  fit <- .C("GLMM_longitDA2",
            nonSilent      = as.integer(!silent),
            Y_c            = as.double(Cy_c),
            R_c            = as.integer(Rc),
            Y_d            = as.integer(Cy_d),
            R_d            = as.integer(Rd),
            dist           = as.integer(Cdist),
            nClust         = as.integer(nClust),
            I              = as.integer(I),
            n              = as.integer(Cn),
            X              = as.double(CX),
            p              = as.integer(Cp),
            fixedIntcpt    = as.integer(CfixedIntcpt),
            Z              = as.double(CZ),
            q              = as.integer(Cq),
            randIntcpt     = as.integer(CrandomIntcpt),
            shiftScale_b   = as.double(CshiftScale_b),
            distribution_b = as.integer(Cdistribution_b),
            keepMCMC       = as.integer(keepMCMC),
            keep_logf      = as.integer(keep.comp.prob),
            info           = as.integer(info),
            Kmax_b         = as.integer(Kmax_b),
            chsigma_eps    = as.double(chsigma_eps),
            chK_b          = as.integer(chK_b),
            chw_b          = as.double(chw_b),
            chmu_b         = as.double(chmu_b),
            chLi_b         = as.double(chLi_b),
            chalpha        = as.double(chalpha),
            f_marg         = double(I * nClust),
            f_cond         = double(I * nClust),
            f_reff         = double(I * nClust),
            logf_marg      = double(I * nClust),
            logf_cond      = double(I * nClust),
            logf_reff      = double(I * nClust),            
            bpred          = double((sum(Cq) + sum(CrandomIntcpt)) * I),
            logf_marg_i    = double(l_logf_i),
            logf_cond_i    = double(l_logf_i),
            logf_reff_i    = double(l_logf_i),            
            nzero_marg     = integer(I * nClust),            
            nzero_cond     = integer(I * nClust),
            nzero_reff     = integer(I * nClust),            
            err            = as.integer(0))
 #           PACKAGE = thispackage)      

  #### Free some space
  fit$chsigma_eps <- NULL
  fit$chK_b       <- NULL
  fit$chw_b       <- NULL    
  fit$chmu_b      <- NULL
  fit$chLi_b      <- NULL
  fit$chalpha     <- NULL      
  
  #### Object to return
  RET <- list()
  idName <- unique(id_cl[[1]])      
  clName <- paste("Cluster ", 1:nClust, sep = "")
  
  #### Predictions of random effects (by cluster)
  RET$bpred <- list()  
  for (cl in 1:nClust){
    if (dimb_cl[cl] == 0) stop("Programming error, contact AK (no random effects for some cluster).")
    if (cl == 1){
      RET$bpred[[cl]] <- matrix(fit$bpred[1:(dimb_cl[cl] * I)], nrow = I, byrow = TRUE)
    }else{
      cumdb <- cumsum(dimb_cl[1:(cl - 1)])
      RET$bpred[[cl]] <- matrix(fit$bpred[(cumdb * I + 1):((cumdb + dimb_cl[cl]) * I)], nrow = I, byrow = TRUE)
    }    
    rownames(RET$bpred[[cl]]) <- idName
  }
  names(RET$bpred) <- clName
  fit$bpred <- NULL
  
  ### Numbers of zero likelihoods
  RET$nzero_marg <- matrix(fit$nzero_marg, nrow = I, ncol = nClust)
  rownames(RET$nzero_marg) <- idName
  colnames(RET$nzero_marg) <- clName  
  fit$nzero_marg <- NULL
  #
  RET$nzero_cond <- matrix(fit$nzero_cond, nrow = I, ncol = nClust)
  rownames(RET$nzero_cond) <- idName
  colnames(RET$nzero_cond) <- clName  
  fit$nzero_cond <- NULL
  #
  RET$nzero_reff <- matrix(fit$nzero_reff, nrow = I, ncol = nClust)  
  rownames(RET$nzero_reff) <- idName
  colnames(RET$nzero_reff) <- clName
  fit$nzero_reff <- NULL  
  
  ### Input for discrimination
  RET$f_marg <- matrix(fit$f_marg, nrow = I, ncol = nClust)
  rownames(RET$f_marg) <- idName
  colnames(RET$f_marg) <- clName  
  fit$f_marg <- NULL
  #
  RET$f_cond <- matrix(fit$f_cond, nrow = I, ncol = nClust)
  rownames(RET$f_cond) <- idName
  colnames(RET$f_cond) <- clName  
  fit$f_cond <- NULL
  #
  RET$f_reff <- matrix(fit$f_reff, nrow = I, ncol = nClust)
  rownames(RET$f_reff) <- idName
  colnames(RET$f_reff) <- clName  
  fit$f_reff <- NULL  
  
  RET$logf_marg <- matrix(fit$logf_marg, nrow = I, ncol = nClust)
  rownames(RET$logf_marg) <- idName
  colnames(RET$logf_marg) <- clName  
  fit$logf_marg <- NULL
  #
  RET$logf_cond <- matrix(fit$logf_cond, nrow = I, ncol = nClust)
  rownames(RET$logf_cond) <- idName
  colnames(RET$logf_cond) <- clName  
  fit$logf_cond <- NULL
  #
  RET$logf_reff <- matrix(fit$logf_reff, nrow = I, ncol = nClust)
  rownames(RET$logf_reff) <- idName
  colnames(RET$logf_reff) <- clName  
  fit$logf_reff <- NULL  
  
  ### Allocation probabilities
  w.prior.mat <- matrix(rep(w.prior, I), nrow = I, ncol = nClust, byrow = TRUE)
  #
  RET$pi_marg <- w.prior.mat * RET$f_marg
  sum.pi_marg <- matrix(rep(apply(RET$pi_marg, 1, sum), nClust), nrow = I, ncol = nClust)
  RET$pi_marg <- RET$pi_marg / sum.pi_marg
  rownames(RET$pi_marg) <- idName
  colnames(RET$pi_marg) <- clName  
  #
  RET$pi_cond <- w.prior.mat * RET$f_cond
  sum.pi_cond <- matrix(rep(apply(RET$pi_cond, 1, sum), nClust), nrow = I, ncol = nClust)
  RET$pi_cond <- RET$pi_cond / sum.pi_cond
  rownames(RET$pi_cond) <- idName
  colnames(RET$pi_cond) <- clName  
  #
  RET$pi_reff <- w.prior.mat * RET$f_reff
  sum.pi_reff <- matrix(rep(apply(RET$pi_reff, 1, sum), nClust), nrow = I, ncol = nClust)
  RET$pi_reff <- RET$pi_reff / sum.pi_reff
  rownames(RET$pi_reff) <- idName
  colnames(RET$pi_reff) <- clName
  
  ### Kept allocation probabilities and related quantities
  if (keep.comp.prob){
    cumkeepMCMC <- cumsum(keepMCMC)

    ### Lists with one matrix (with sampled f's) for each cluster
    f_marg_i <- f_cond_i <- f_reff_i <- list()
    for (cl in 1:nClust){
      if (cl == 1){        
        f_marg_i[[cl]] <- matrix(exp(fit$logf_marg_i[1:(I * keepMCMC[1])]), ncol = I, nrow = keepMCMC[1], byrow = TRUE)
        colnames(f_marg_i[[cl]]) <- idName
        #
        f_cond_i[[cl]] <- matrix(exp(fit$logf_cond_i[1:(I * keepMCMC[1])]), ncol = I, nrow = keepMCMC[1], byrow = TRUE)
        colnames(f_cond_i[[cl]]) <- idName
        #
        f_reff_i[[cl]] <- matrix(exp(fit$logf_reff_i[1:(I * keepMCMC[1])]), ncol = I, nrow = keepMCMC[1], byrow = TRUE)
        colnames(f_reff_i[[cl]]) <- idName                
      }else{
        f_marg_i[[cl]] <- matrix(exp(fit$logf_marg_i[(I * cumkeepMCMC[cl - 1] + 1):(I * cumkeepMCMC[cl])]), ncol = I, nrow = keepMCMC[cl], byrow = TRUE)
        colnames(f_marg_i[[cl]]) <- idName
        #
        f_cond_i[[cl]] <- matrix(exp(fit$logf_cond_i[(I * cumkeepMCMC[cl - 1] + 1):(I * cumkeepMCMC[cl])]), ncol = I, nrow = keepMCMC[cl], byrow = TRUE)
        colnames(f_cond_i[[cl]]) <- idName
        #
        f_reff_i[[cl]] <- matrix(exp(fit$logf_reff_i[(I * cumkeepMCMC[cl - 1] + 1):(I * cumkeepMCMC[cl])]), ncol = I, nrow = keepMCMC[cl], byrow = TRUE)
        colnames(f_reff_i[[cl]]) <- idName        
      }    
    }
    fit$logf_marg_i <- NULL
    fit$logf_cond_i <- NULL
    fit$logf_reff_i <- NULL

    ### Put everything into matrices
    ### columns: i = 1, cluster 1, 2, ..., nClust; i = 2, cluster 1, 2, ..., nClust, etc.
    ### rows:    MCMC iterations (if different lengths, take first maxkeepMCMC iterations)
    maxkeepMCMC <- max(keepMCMC)
    RET$pi_marg_i <- RET$pi_cond_i <- RET$pi_reff_i <- matrix(NA, ncol = I * nClust, nrow = maxkeepMCMC)
    colnames(RET$pi_marg_i) <- colnames(RET$pi_cond_i) <- colnames(RET$pi_reff_i) <- paste(rep(idName, each = nClust), ": ", rep(1:nClust, I), sep = "")
    #
    for (cl in 1:nClust){
      f_marg_i[[cl]] <- f_marg_i[[cl]] * w.prior[cl]
      f_cond_i[[cl]] <- f_cond_i[[cl]] * w.prior[cl]
      f_reff_i[[cl]] <- f_reff_i[[cl]] * w.prior[cl]      
    }    
    
    for (i in 1:I){
      for (cl in 1:nClust){
        RET$pi_marg_i[, (i - 1)*nClust + cl] <- f_marg_i[[cl]][1:maxkeepMCMC, i]        
        RET$pi_cond_i[, (i - 1)*nClust + cl] <- f_cond_i[[cl]][1:maxkeepMCMC, i]
        RET$pi_reff_i[, (i - 1)*nClust + cl] <- f_reff_i[[cl]][1:maxkeepMCMC, i]        
      }
      sum_pi_marg_i <- apply(RET$pi_marg_i[, ((i - 1)*nClust + 1):(i*nClust)], 1, sum, na.rm = FALSE)
      RET$pi_marg_i[, ((i - 1)*nClust + 1):(i*nClust)] <- RET$pi_marg_i[, ((i - 1)*nClust + 1):(i*nClust)] / matrix(rep(sum_pi_marg_i, nClust), nrow = maxkeepMCMC, ncol = nClust, byrow = FALSE)
      #
      sum_pi_cond_i <- apply(RET$pi_cond_i[, ((i - 1)*nClust + 1):(i*nClust)], 1, sum, na.rm = FALSE)
      RET$pi_cond_i[, ((i - 1)*nClust + 1):(i*nClust)] <- RET$pi_cond_i[, ((i - 1)*nClust + 1):(i*nClust)] / matrix(rep(sum_pi_cond_i, nClust), nrow = maxkeepMCMC, ncol = nClust, byrow = FALSE)
      #
      sum_pi_reff_i <- apply(RET$pi_reff_i[, ((i - 1)*nClust + 1):(i*nClust)], 1, sum, na.rm = FALSE)
      RET$pi_reff_i[, ((i - 1)*nClust + 1):(i*nClust)] <- RET$pi_reff_i[, ((i - 1)*nClust + 1):(i*nClust)] / matrix(rep(sum_pi_reff_i, nClust), nrow = maxkeepMCMC, ncol = nClust, byrow = FALSE)      
    }
    rm(list = c("f_marg_i", "f_cond_i", "f_reff_i"))

    ### Calculate posterior means, medians and HPD's
    RET$pi_marg_mean <- matrix(apply(RET$pi_marg_i, 2, mean, na.rm = TRUE), ncol = nClust, nrow = I, byrow = TRUE)
    RET$pi_marg_med  <- matrix(apply(RET$pi_marg_i, 2, median, na.rm = TRUE), ncol = nClust, nrow = I, byrow = TRUE)
    hpd_marg <- coda::HPDinterval(coda::mcmc(RET$pi_marg_i), prob = level)
    RET$pi_marg_low <- matrix(hpd_marg[, "lower"], ncol = nClust, nrow = I, byrow = TRUE)
    RET$pi_marg_upp <- matrix(hpd_marg[, "upper"], ncol = nClust, nrow = I, byrow = TRUE)
    rownames(RET$pi_marg_mean) <- rownames(RET$pi_marg_med) <- rownames(RET$pi_marg_low) <- rownames(RET$pi_marg_upp) <- idName
    colnames(RET$pi_marg_mean) <- colnames(RET$pi_marg_med) <- colnames(RET$pi_marg_low) <- colnames(RET$pi_marg_upp) <- clName
    #     
    RET$pi_cond_mean <- matrix(apply(RET$pi_cond_i, 2, mean, na.rm = TRUE), ncol = nClust, nrow = I, byrow = TRUE)
    RET$pi_cond_med  <- matrix(apply(RET$pi_cond_i, 2, median, na.rm = TRUE), ncol = nClust, nrow = I, byrow = TRUE)
    hpd_cond <- coda::HPDinterval(coda::mcmc(RET$pi_cond_i), prob = level)
    RET$pi_cond_low <- matrix(hpd_cond[, "lower"], ncol = nClust, nrow = I, byrow = TRUE)
    RET$pi_cond_upp <- matrix(hpd_cond[, "upper"], ncol = nClust, nrow = I, byrow = TRUE)
    rownames(RET$pi_cond_mean) <- rownames(RET$pi_cond_med) <- rownames(RET$pi_cond_low) <- rownames(RET$pi_cond_upp) <- idName
    colnames(RET$pi_cond_mean) <- colnames(RET$pi_cond_med) <- colnames(RET$pi_cond_low) <- colnames(RET$pi_cond_upp) <- clName
    #
    RET$pi_reff_mean <- matrix(apply(RET$pi_reff_i, 2, mean, na.rm = TRUE), ncol = nClust, nrow = I, byrow = TRUE)
    RET$pi_reff_med  <- matrix(apply(RET$pi_reff_i, 2, median, na.rm = TRUE), ncol = nClust, nrow = I, byrow = TRUE)
    hpd_reff <- coda::HPDinterval(coda::mcmc(RET$pi_reff_i), prob = level)
    RET$pi_reff_low <- matrix(hpd_reff[, "lower"], ncol = nClust, nrow = I, byrow = TRUE)
    RET$pi_reff_upp <- matrix(hpd_reff[, "upper"], ncol = nClust, nrow = I, byrow = TRUE)
    rownames(RET$pi_reff_mean) <- rownames(RET$pi_reff_med) <- rownames(RET$pi_reff_low) <- rownames(RET$pi_reff_upp) <- idName
    colnames(RET$pi_reff_mean) <- colnames(RET$pi_reff_med) <- colnames(RET$pi_reff_low) <- colnames(RET$pi_reff_upp) <- clName        
  }  
  
  return(RET)
}

