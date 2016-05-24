##
##  PURPOSE:   Longitudinal discriminant analysis
##             based on GLMM MCMC fits (with possibly several response variables)
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  LOG:       05/08/2009: created as GLMM_longitClust, only allows for continuous responses
##             28/10/2009: renamed to GLMM_longitDA, still allows for continuous responses only
##             16/04/2015: minor revision (stop added when some response variable is not continuous)
##
##  FUNCTIONS:  GLMM_longitDA
##
## ================================================================================================

## *************************************************************
## GLMM_longitDA
## *************************************************************
##
GLMM_longitDA <- function(mod, w.prior, y, id, time, x, z, xz.common = TRUE, info)
{
  thispackage <- "mixAK" 
  silent <- TRUE
  
########## ========== Some input checks ========== ##########
########## ======================================= ##########
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
  if (!xz.common){
    if (!is.list(x)) stop("x must be a list when xz.common is FALSE")
    if (length(x) != nClust) stop(paste("x must be a list of length ", nClust, sep=""))
    if (!is.list(z)) stop("z must be a list when xz.common is FALSE")
    if (length(z) != nClust) stop(paste("z must be a list of length ", nClust, sep=""))    
  }  


########## ========== Work out the data ============ ##########
########## ========================================= ##########
  Cp <- Cq <- CfixedIntcpt <- CrandomIntcpt <- numeric()
  CX <- CZ <- numeric()
  ## CXtX <- CZitZi <- numeric()                              ### CODE REMOVED ON 02/11/2009
  CshiftScale_b <- numeric()
  Kmax_b <- numeric()
  chK_b <- chw_b <- chmu_b <- chLi_b <- numeric()
  chsigma_eps <- numeric()
  chalpha <- numeric()
  
  for (cl in 1:nClust){     ### Loop over clusters

    ##### Design stuff
    ##### ---------------
    if (xz.common){
      dd <- GLMM_MCMCdata(y=y, dist=mod[[cl]]$dist, id=id, time=time, x=x, z=z, random.intercept=mod[[cl]]$random.intercept)
    }else{
      dd <- GLMM_MCMCdata(y=y, dist=mod[[cl]]$dist, id=id, time=time, x=x[[cl]], z=z[[cl]], random.intercept=mod[[cl]]$random.intercept)
    }

    ifit <- GLMM_MCMCifit(do.init=FALSE, na.complete=TRUE,
                          y=dd$y, dist=dd$dist, id=dd$id, time=dd$time, x=dd$x, z=dd$z, random.intercept=dd$random.intercept,
                          xempty=dd$xempty, zempty=dd$zempty, Rc=dd$Rc, Rd=dd$Rd,
                          p=dd$p, p_fi=dd$p_fi, q=dd$q, q_ri=dd$q_ri, lalpha=dd$lalpha, dimb=dd$dimb)
    ## ifit$x[[s]], ifit$z[[s]] now contains also intercept columns if these should be included
    ## * also rows corresponding to missing values are removed
    ## * there is the same number of observations for each cluster
    
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
      TIME <- ifit$time

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
    
    ##### Check whether supplied x and z is compatible with these used to fit mod
    ##### Double check that there is the same number of observations for each response
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

      if (any(ifit$n[[s]] != ifit$n[[1]])) stop("BUG in the function, contact AK")
    }

    ### ===== CODE REMOVED ON 02/11/2009 =====
    ##### Calculate CZitZi to be used here (it is different from that needed by GLMM_MCMC function)
    ##### and now included in ifit$CZitZi
    ##### REMARK:  Due to a loop, this is relatively slow for larger data, do it in C++ in a future???
    ##ZitZi <- numeric(0)
    ##cumn <- c(0, cumsum(ifit$n[[1]]))
    ##for (i in 1:I){                       ## loop over longitudinal profiles
    ##  for (j in 1:ifit$n[[1]][i]){        ## loop over observations within longitudinal profiles
    ##    for (s in 1:R){
    ##      if (dd$q_ri[s]){
    ##        tmpZ <- matrix(ifit$z[[s]][(cumn[i]+1):(cumn[i]+j),], ncol=ncol(ifit$z[[s]]))
    ##        tmpZtZ <- t(tmpZ) %*% tmpZ
    ##        ZitZi <- c(ZitZi, tmpZtZ[lower.tri(tmpZtZ, diag=TRUE)])
    ##      }  
    ##    }          
    ##  }  
    ##}
    ##cat("length ZitZi[", cl, "] = ", length(ZitZi), "\n", sep="")
    ### ===== END OF CODE REMOVED ON 02/11/2009 =====
    
    ##### Vectors to be supplied to C++
    ##### ---------------------------------
    Cp <- c(Cp, dd$p)
    Cq <- c(Cq, dd$q)    
    CfixedIntcpt  <- c(CfixedIntcpt, dd$CfixedIntcpt)
    CrandomIntcpt <- c(CrandomIntcpt, dd$CrandomIntcpt)
    CX     <- c(CX, ifit$CX)
    ##CXtX   <- c(CXtX, ifit$CXtX)           ## CODE REMOVED ON 02/11/2009
    CZ     <- c(CZ, ifit$CZ)
    ##CZitZi <- c(CZitZi, ZitZi)             ## CODE REMOVED ON 02/11/2009

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
  
  ### !!! CX, CXtX, CZ, CZitZi contain one 0 when there is no such matrix
  ###     in a model for specific cluster

  if (missing(info)) info <- min(keepMCMC)
  if (info <= 0) info <- 1

  ##### ===============================================================
  #####                                                           #####
  ##### Validation of C++ code of GLMM_longitDA function          #####
  ##### may follow here.                                          #####
  ##### Corresponding code is available in                        #####
  ##### ~/Rlib/mixAK/Rtemp/GLMM_longitDA/GLMM_longitDA_validC.R   #####
  #####                                                           #####
  ##### ===============================================================

  if (FALSE){
    cat("Cy_c:\n")
    print(Cy_c)

    cat("Cy_d:\n")
    print(Cy_d)

    cat("CX:\n")
    print(CX)

    cat("CZ:\n")
    print(CZ)

    ##cat("CZitZi:\n")
    ##print(CZitZi)

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
  
  if (any(Cdist != 0)){
      
    stop("Not (yet) implemented if any of the response variables is not continuous.")
    
  }else{
      
    fit <- .C("GLMM_longitDA",
              Y_c          = as.double(Cy_c),
              R_c          = as.integer(Rc),
              Y_d          = as.integer(Cy_d),
              R_d          = as.integer(Rd),
              dist         = as.integer(Cdist),
              nClust       = as.integer(nClust),
              I            = as.integer(I),
              n            = as.integer(Cn),
              X            = as.double(CX),
              p            = as.integer(Cp),
              fixedIntcpt  = as.integer(CfixedIntcpt),
              Z            = as.double(CZ),
              #SZitZiS      = as.double(CZitZi),                  ### REMOVED ON 02/11/2009
              q            = as.integer(Cq),
              randIntcpt   = as.integer(CrandomIntcpt),
              shiftScale_b = as.double(CshiftScale_b),
              keepMCMC     = as.integer(keepMCMC),
              info         = as.integer(info),
              Kmax_b       = as.integer(Kmax_b),
              chsigma_eps  = as.double(chsigma_eps),
              chK_b        = as.integer(chK_b),
              chw_b        = as.double(chw_b),
              chmu_b       = as.double(chmu_b),
              chLi_b       = as.double(chLi_b),
              chalpha      = as.double(chalpha),
              pi_marg      = double(nrow_pi * nClust),
              pi_cond      = double(nrow_pi * nClust),
              pi_ranef     = double(nrow_pi * nClust),
              err          = as.integer(0),
              PACKAGE = thispackage)
  }  
  
  if (nrow_pi == 1){
    fit$pi_ranef <- w.prior * fit$pi_ranef
    fit$pi_ranef <- fit$pi_ranef / sum(fit$pi_ranef)

    fit$pi_cond <- w.prior * fit$pi_cond
    fit$pi_cond <- fit$pi_cond / sum(fit$pi_cond)

    fit$pi_marg <- w.prior * fit$pi_marg
    fit$pi_marg <- fit$pi_marg / sum(fit$pi_marg)    
  }else{     
    w.priorMat <- matrix(rep(w.prior, nrow_pi), ncol = nClust, byrow = TRUE)
    #
    fit$pi_ranef <- w.priorMat * matrix(fit$pi_ranef, ncol = nClust)
    fit$pi_ranef <- fit$pi_ranef / apply(fit$pi_ranef, 1, sum)
    #
    fit$pi_cond <- w.priorMat * matrix(fit$pi_cond, ncol = nClust)
    fit$pi_cond <- fit$pi_cond / apply(fit$pi_cond, 1, sum)
    #
    fit$pi_marg <- w.priorMat * matrix(fit$pi_marg, ncol = nClust)
    fit$pi_marg <- fit$pi_marg / apply(fit$pi_marg, 1, sum)
  }  

  RET <- list(ident = data.frame(id = ID, time = TIME),
              marg  = fit$pi_marg,
              cond  = fit$pi_cond,
              ranef = fit$pi_ranef)
  
  return(RET)  
}
